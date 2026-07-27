require 'json'
require 'set'

module AGS

  #{{{ AGNOSTIC GRN HELPERS

  helper :grn_all_contexts do
    contexts = []
    TREATMENTS.each do |treatment|
      TIME_POINTS.each do |time_point|
        contexts << [treatment, time_point]
      end
    end
    contexts
  end

  helper :grn_context_name do |treatment, time_point|
    [treatment, time_point] * "-"
  end

  helper :grn_sign_value do |value, threshold = 0.0|
    return 0 if value.nil?
    value = value.to_f
    return 1 if value > threshold
    return -1 if value < -threshold
    0
  end

  helper :grn_state_label do |sign|
    case sign.to_i
    when 1 then "positive"
    when -1 then "negative"
    else "absent"
    end
  end

  helper :grn_entity_states do |grn, entity_type, sign_threshold = 0.0|
    states = {}
    case entity_type.to_s
    when "tg", "target", "target_state"
      grn[:tgs].each do |gene, direction|
        states[gene] = direction == "increase" ? 1 : -1
      end
    when "tf", "tf_activity"
      grn[:tfs].each do |tf, activity|
        states[tf] = grn_sign_value(activity, sign_threshold)
      end
    when "target_pressure"
      pressure = Hash.new(0.0)
      grn[:edges].each do |edge_key, sign|
        tf, tg = edge_key.split("~")
        next unless tf && tg
        next unless grn[:tfs].key?(tf)
        pressure[tg] += grn[:tfs][tf].to_f * sign.to_i
      end
      pressure.each do |tg, value|
        states[tg] = grn_sign_value(value, sign_threshold)
      end
    when "edge", "edge_pressure"
      grn[:edges].each do |edge_key, sign|
        tf, _tg = edge_key.split("~")
        next unless tf && grn[:tfs].key?(tf)
        states[edge_key] = grn_sign_value(grn[:tfs][tf].to_f * sign.to_i, sign_threshold)
      end
    else
      raise "Unknown GRN entity_type: #{entity_type}"
    end
    states
  end

  helper :grn_entity_values do |grn, matrix_type, value_mode = "sign", sign_threshold = 0.0|
    values = {}
    case matrix_type.to_s
    when "tg", "target", "target_state"
      grn[:tgs].each do |gene, direction|
        sign = direction == "increase" ? 1 : -1
        values[gene] = value_mode == "label" ? grn_state_label(sign) : sign
      end
    when "tf", "tf_activity"
      grn[:tfs].each do |tf, activity|
        if value_mode == "value"
          values[tf] = activity.to_f
        else
          sign = grn_sign_value(activity, sign_threshold)
          values[tf] = value_mode == "label" ? grn_state_label(sign) : sign
        end
      end
    when "target_pressure"
      pressure = Hash.new(0.0)
      grn[:edges].each do |edge_key, sign|
        tf, tg = edge_key.split("~")
        next unless tf && tg
        next unless grn[:tfs].key?(tf)
        pressure[tg] += grn[:tfs][tf].to_f * sign.to_i
      end
      pressure.each do |tg, val|
        if value_mode == "value"
          values[tg] = val.round(4)
        else
          sign = grn_sign_value(val, sign_threshold)
          values[tg] = value_mode == "label" ? grn_state_label(sign) : sign
        end
      end
    when "edge", "edge_pressure"
      grn[:edges].each do |edge_key, sign|
        tf, _tg = edge_key.split("~")
        next unless tf && grn[:tfs].key?(tf)
        val = grn[:tfs][tf].to_f * sign.to_i
        if value_mode == "value"
          values[edge_key] = val.round(4)
        else
          s = grn_sign_value(val, sign_threshold)
          values[edge_key] = value_mode == "label" ? grn_state_label(s) : s
        end
      end
    else
      raise "Unknown GRN matrix_type: #{matrix_type}"
    end
    values
  end

  helper :grn_transition_label do |state1, state2|
    [grn_state_label(state1), grn_state_label(state2)] * "-to-"
  end

  helper :grn_neutral_component_class do |combo_state, c1_state, c2_state|
    combo_present = combo_state.to_i != 0
    c1_present = c1_state.to_i != 0
    c2_present = c2_state.to_i != 0

    return "absent_all" unless combo_present || c1_present || c2_present
    return "combo_specific_present" if combo_present && !c1_present && !c2_present
    return "component_specific_absent_in_combo" if !combo_present && (c1_present || c2_present)

    if combo_present
      same1 = c1_present && combo_state.to_i == c1_state.to_i
      same2 = c2_present && combo_state.to_i == c2_state.to_i
      inv1 = c1_present && combo_state.to_i == -c1_state.to_i
      inv2 = c2_present && combo_state.to_i == -c2_state.to_i
      return "shared_same" if same1 && same2
      return "component1_like" if same1 && !same2
      return "component2_like" if same2 && !same1
      return "combo_sign_inverted" if inv1 || inv2
      return "mixed_present"
    end

    "mixed_present"
  end

  #{{{ AGNOSTIC STATE MATRIX

  dep :grns, jobname: 'Default'
  input :matrix_type, :select, "Matrix type", "target_state", :select_options => %w(target_state tf_activity target_pressure edge_pressure), :required => true
  input :value_mode, :select, "Value mode", "sign", :select_options => %w(sign value label), :required => true
  input :sign_threshold, :float, "Threshold for sign calls on activity/pressure", 0.0
  task :grn_state_matrix => :tsv do |matrix_type, value_mode, sign_threshold|
    contexts = grn_all_contexts
    fields = contexts.collect { |treatment, time_point| grn_context_name(treatment, time_point) }
    rows = Hash.new { |h, k| h[k] = {} }

    contexts.each do |treatment, time_point|
      context = grn_context_name(treatment, time_point)
      grn = grn_load_json(treatment, time_point)
      values = grn_entity_values(grn, matrix_type, value_mode, sign_threshold)
      values.each do |entity, value|
        rows[entity][context] = value
      end
    end

    result = TSV.setup({}, :key_field => "Entity", :fields => fields, :type => :list)
    rows.keys.sort.each do |entity|
      result[entity] = fields.collect { |field| rows[entity].fetch(field, value_mode == "label" ? "absent" : 0).to_s }
    end
    result
  end

  #{{{ AGNOSTIC SIGNATURE MODULES

  dep :grns, jobname: 'Default'
  input :matrix_type, :select, "Matrix type", "target_state", :select_options => %w(target_state tf_activity target_pressure edge_pressure), :required => true
  input :min_size, :integer, "Minimum number of entities per discovered signature", 5
  input :sign_threshold, :float, "Threshold for sign calls on activity/pressure", 0.0
  task :grn_signature_modules => :tsv do |matrix_type, min_size, sign_threshold|
    min_size = min_size.to_i
    min_size = 5 if min_size <= 0
    contexts = grn_all_contexts
    context_names = contexts.collect { |treatment, time_point| grn_context_name(treatment, time_point) }
    rows = Hash.new { |h, k| h[k] = Array.new(contexts.length, 0) }

    contexts.each_with_index do |(treatment, time_point), i|
      grn = grn_load_json(treatment, time_point)
      states = grn_entity_states(grn, matrix_type, sign_threshold)
      states.each do |entity, state|
        rows[entity][i] = state.to_i
      end
    end

    groups = Hash.new { |h, k| h[k] = [] }
    rows.each do |entity, signature|
      next if signature.all? { |v| v.to_i == 0 }
      groups[signature * ","] << entity
    end

    fields = %w(Size Positive_contexts Negative_contexts Nonzero_contexts Signature Entities)
    result = TSV.setup({}, :key_field => "Discovered_module", :fields => fields, :type => :list)
    mod_i = 0
    groups.sort_by { |_sig, entities| [-entities.length, _sig] }.each do |sig, entities|
      next if entities.length < min_size
      signature = sig.split(',').collect(&:to_i)
      mod_i += 1
      pos = []
      neg = []
      nonzero = []
      signature.each_with_index do |v, i|
        if v > 0
          pos << context_names[i]
          nonzero << context_names[i]
        elsif v < 0
          neg << context_names[i]
          nonzero << context_names[i]
        end
      end
      result["M%03d" % mod_i] = [
        entities.length.to_s,
        pos * ",",
        neg * ",",
        nonzero * ",",
        sig,
        grn_format_gene_list(entities, 200)
      ]
    end
    result
  end

  #{{{ AGNOSTIC COMBINATION CONTRAST

  dep :grns, jobname: 'Default'
  input :combo_treatment, :select, "Combination treatment", "INT_PD_PI", :select_options => TREATMENTS, :required => true
  input :component1, :select, "Component 1 treatment", "PI", :select_options => TREATMENTS, :required => true
  input :component2, :select, "Component 2 treatment", "PD", :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", 8, :select_options => TIME_POINTS, :required => true
  input :entity_type, :select, "Entity type", "target_state", :select_options => %w(target_state tf_activity target_pressure edge_pressure), :required => true
  input :sign_threshold, :float, "Threshold for sign calls on activity/pressure", 0.0
  task :grn_neutral_contrast => :tsv do |combo_treatment, component1, component2, time_point, entity_type, sign_threshold|
    combo = grn_entity_states(grn_load_json(combo_treatment, time_point), entity_type, sign_threshold)
    c1 = grn_entity_states(grn_load_json(component1, time_point), entity_type, sign_threshold)
    c2 = grn_entity_states(grn_load_json(component2, time_point), entity_type, sign_threshold)
    entities = Set.new(combo.keys) | Set.new(c1.keys) | Set.new(c2.keys)

    fields = %w(Combo_state Component1_state Component2_state Combo_label Component1_label Component2_label Classification)
    result = TSV.setup({}, :key_field => "Entity", :fields => fields, :type => :list)
    entities.sort.each do |entity|
      cs = combo[entity].to_i
      s1 = c1[entity].to_i
      s2 = c2[entity].to_i
      result[entity] = [
        cs.to_s, s1.to_s, s2.to_s,
        grn_state_label(cs), grn_state_label(s1), grn_state_label(s2),
        grn_neutral_component_class(cs, s1, s2)
      ]
    end
    result
  end

  #{{{ AGNOSTIC TRANSITION CONTRAST

  dep :grns, jobname: 'Default'
  input :combo_treatment, :select, "Combination treatment", "INT_PD_PI", :select_options => TREATMENTS, :required => true
  input :component1, :select, "Component 1 treatment", "PI", :select_options => TREATMENTS, :required => true
  input :component2, :select, "Component 2 treatment", "PD", :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Earlier timepoint", 4, :select_options => TIME_POINTS, :required => true
  input :time_point2, :select, "Later timepoint", 8, :select_options => TIME_POINTS, :required => true
  input :entity_type, :select, "Entity type", "target_state", :select_options => %w(target_state tf_activity target_pressure edge_pressure), :required => true
  input :sign_threshold, :float, "Threshold for sign calls on activity/pressure", 0.0
  task :grn_neutral_transition_contrast => :tsv do |combo_treatment, component1, component2, time_point, time_point2, entity_type, sign_threshold|
    c_a = grn_entity_states(grn_load_json(combo_treatment, time_point), entity_type, sign_threshold)
    c_b = grn_entity_states(grn_load_json(combo_treatment, time_point2), entity_type, sign_threshold)
    c1_a = grn_entity_states(grn_load_json(component1, time_point), entity_type, sign_threshold)
    c1_b = grn_entity_states(grn_load_json(component1, time_point2), entity_type, sign_threshold)
    c2_a = grn_entity_states(grn_load_json(component2, time_point), entity_type, sign_threshold)
    c2_b = grn_entity_states(grn_load_json(component2, time_point2), entity_type, sign_threshold)
    entities = Set.new(c_a.keys) | Set.new(c_b.keys) | Set.new(c1_a.keys) | Set.new(c1_b.keys) | Set.new(c2_a.keys) | Set.new(c2_b.keys)

    fields = %w(Combo_transition Component1_transition Component2_transition Combo_T1 Combo_T2 Component1_T2 Component2_T2 Classification)
    result = TSV.setup({}, :key_field => "Entity", :fields => fields, :type => :list)
    entities.sort.each do |entity|
      combo_transition = grn_transition_label(c_a[entity].to_i, c_b[entity].to_i)
      c1_transition = grn_transition_label(c1_a[entity].to_i, c1_b[entity].to_i)
      c2_transition = grn_transition_label(c2_a[entity].to_i, c2_b[entity].to_i)
      classification = if combo_transition == c1_transition && combo_transition == c2_transition
                         "shared_transition"
                       elsif combo_transition != c1_transition && combo_transition != c2_transition
                         "combo_distinct_transition"
                       elsif combo_transition == c1_transition
                         "component1_like_transition"
                       elsif combo_transition == c2_transition
                         "component2_like_transition"
                       else
                         "mixed_transition"
                       end
      result[entity] = [
        combo_transition,
        c1_transition,
        c2_transition,
        grn_state_label(c_a[entity].to_i),
        grn_state_label(c_b[entity].to_i),
        grn_state_label(c1_b[entity].to_i),
        grn_state_label(c2_b[entity].to_i),
        classification
      ]
    end
    result
  end


  #{{{ AGNOSTIC GLOBAL MOTIFS

  dep :grns, jobname: 'Default'
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  input :max_motifs, :integer, "Maximum motifs to return", 500
  input :only_coherent, :boolean, "Only return motifs where both TFs coherently explain the target", true
  task :grn_global_motifs => :tsv do |treatment, time_point, max_motifs, only_coherent|
    grn = grn_load_json(treatment, time_point)
    max_motifs = max_motifs.to_i
    max_motifs = 500 if max_motifs <= 0

    tfs = grn[:tfs]
    tgs = grn[:tgs]
    edges = grn[:edges]

    outgoing = Hash.new { |h, k| h[k] = {} }
    edges.each do |edge_key, sign|
      tf, tg = edge_key.split("~")
      next unless tf && tg
      next unless tfs.key?(tf)
      outgoing[tf][tg] = sign.to_i
    end

    fields = %w(
      Motif_class TF_A TF_B Target A_activity B_activity
      A_to_B_sign A_to_Target_sign B_to_Target_sign Target_direction B_direction
      A_to_B_pressure A_to_Target_pressure B_to_Target_pressure
      Target_coherence B_transcript_coherence Motif_score Interpretation_hint
    )
    result = TSV.setup({}, :key_field => "Motif", :fields => fields, :type => :list)

    motifs = []
    outgoing.each do |tf_a, targets_a|
      a_activity = tfs[tf_a].to_f
      targets_a.each do |tf_b, sign_ab|
        next if tf_b == tf_a
        next unless tfs.key?(tf_b)
        next unless outgoing.key?(tf_b)
        shared_targets = targets_a.keys & outgoing[tf_b].keys & tgs.keys
        shared_targets.each do |target|
          next if target == tf_a || target == tf_b
          sign_at = targets_a[target]
          sign_bt = outgoing[tf_b][target]
          b_activity = tfs[tf_b].to_f
          p_ab = a_activity * sign_ab
          p_at = a_activity * sign_at
          p_bt = b_activity * sign_bt
          target_direction = tgs[target]
          b_direction = tgs[tf_b]

          target_coherence_count = 0
          target_coherence_count += 1 if grn_coherent?(p_at, target_direction)
          target_coherence_count += 1 if grn_coherent?(p_bt, target_direction)
          next if only_coherent && target_coherence_count < 2

          target_coherence = case target_coherence_count
                             when 2 then "both_coherent"
                             when 1 then "one_coherent"
                             else "conflict"
                             end
          b_transcript_coherence = if b_direction.nil?
                                     "b_not_dynamic_tg"
                                   elsif grn_coherent?(p_ab, b_direction)
                                     "coherent"
                                   else
                                     "conflict"
                                   end

          motif_class = if p_at > 0 && p_bt > 0
                          "coherent_activation_ffl"
                        elsif p_at < 0 && p_bt < 0
                          "coherent_repression_ffl"
                        else
                          "incoherent_ffl"
                        end

          motif_score = p_at.abs + p_bt.abs + (target_coherence_count * 2.0)
          motif_score += p_ab.abs * 0.25
          motif_score += 1.0 if b_transcript_coherence == "coherent"
          hint = if motif_class == "coherent_repression_ffl" && target_direction == "decrease"
                   "Global feed-forward repression coherently explains target decrease."
                 elsif motif_class == "coherent_activation_ffl" && target_direction == "increase"
                   "Global feed-forward activation coherently explains target increase."
                 elsif motif_class == "incoherent_ffl"
                   "Global incoherent feed-forward motif may buffer or create delayed reversal."
                 else
                   "Global feed-forward motif has mixed relation to observed target direction."
                 end

          motifs << [motif_score, tf_a, tf_b, target, [
            motif_class,
            tf_a,
            tf_b,
            target,
            a_activity.round(4).to_s,
            b_activity.round(4).to_s,
            sign_ab.to_s,
            sign_at.to_s,
            sign_bt.to_s,
            target_direction.to_s,
            b_direction.nil? ? "absent" : b_direction.to_s,
            p_ab.round(4).to_s,
            p_at.round(4).to_s,
            p_bt.round(4).to_s,
            target_coherence,
            b_transcript_coherence,
            motif_score.round(4).to_s,
            hint
          ]]
        end
      end
    end

    motifs.sort_by { |score, _a, _b, _t, _values| -score }.first(max_motifs).each_with_index do |(_score, tf_a, tf_b, target, values), i|
      result[[i + 1, tf_a, tf_b, target] * "~"] = values
    end
    result
  end


  #{{{ AGNOSTIC SIGNATURE ANNOTATION

  dep :grns, jobname: 'Default'
  input :matrix_type, :select, "Matrix type", "target_state", :select_options => %w(target_state target_pressure tf_activity edge_pressure), :required => true
  input :min_size, :integer, "Minimum number of entities per discovered signature", 5
  input :sign_threshold, :float, "Threshold for sign calls on activity/pressure", 0.0
  input :top_n, :integer, "Number of top annotations to report", 10
  task :grn_annotate_signatures => :tsv do |matrix_type, min_size, sign_threshold, top_n|
    min_size = min_size.to_i
    min_size = 5 if min_size <= 0
    top_n = top_n.to_i
    top_n = 10 if top_n <= 0

    contexts = grn_all_contexts
    context_names = contexts.collect { |treatment, time_point| grn_context_name(treatment, time_point) }
    curated_modules = grn_load_modules
    curated_gene_modules = grn_gene_modules(curated_modules)

    rows = Hash.new { |h, k| h[k] = Array.new(contexts.length, 0) }
    grns_by_context = {}

    contexts.each_with_index do |(treatment, time_point), i|
      grn = grn_load_json(treatment, time_point)
      grns_by_context[context_names[i]] = grn
      states = grn_entity_states(grn, matrix_type, sign_threshold)
      states.each do |entity, state|
        rows[entity][i] = state.to_i
      end
    end

    groups = Hash.new { |h, k| h[k] = [] }
    rows.each do |entity, signature|
      next if signature.all? { |v| v.to_i == 0 }
      groups[signature * ","] << entity
    end

    fields = %w(
      Size Positive_contexts Negative_contexts Nonzero_contexts
      Curated_module_overlap Top_coherent_TFs Top_positive_TFs Top_negative_TFs
      Top_TF_contexts Top_TF_targets Representative_entities Signature
    )
    result = TSV.setup({}, :key_field => "Discovered_module", :fields => fields, :type => :list)

    mod_i = 0
    groups.sort_by { |_sig, entities| [-entities.length, _sig] }.each do |sig, entities|
      next if entities.length < min_size
      signature = sig.split(',').collect(&:to_i)
      mod_i += 1
      pos_contexts = []
      neg_contexts = []
      nonzero_contexts = []
      signature.each_with_index do |v, i|
        if v > 0
          pos_contexts << context_names[i]
          nonzero_contexts << context_names[i]
        elsif v < 0
          neg_contexts << context_names[i]
          nonzero_contexts << context_names[i]
        end
      end

      # Post-hoc overlap with curated modules; this is annotation only, not used
      # to define the discovered groups.
      overlap_counts = Hash.new(0)
      entities.each do |entity|
        gene = entity.to_s.split('~').last
        (curated_gene_modules[gene] || []).each { |mod| overlap_counts[mod] += 1 }
      end
      curated_overlap = overlap_counts.sort_by { |mod, count| [-count, mod] }.first(top_n).collect { |mod, count| "#{mod}:#{count}" } * ","

      # For gene-like signatures, collect coherent TF pressure after discovery.
      tf_abs = Hash.new(0.0)
      tf_pos = Hash.new(0.0)
      tf_neg = Hash.new(0.0)
      tf_contexts = Hash.new { |h, k| h[k] = Set.new }
      tf_targets = Hash.new { |h, k| h[k] = Set.new }
      entity_set = Set.new(entities.collect { |e| e.to_s })

      nonzero_contexts.each do |context|
        grn = grns_by_context[context]
        grn[:edges].each do |edge_key, edge_sign|
          tf, tg = edge_key.split("~")
          next unless tf && tg
          # For edge signatures, annotate by edge source/target; for all other
          # signatures, use entity as target gene/TF name.
          in_group = if matrix_type.to_s == "edge_pressure"
                       entity_set.include?(edge_key)
                     else
                       entity_set.include?(tg)
                     end
          next unless in_group
          next unless grn[:tfs].key?(tf)
          tg_direction = grn[:tgs][tg]
          pressure = grn[:tfs][tf].to_f * edge_sign.to_i
          next unless grn_coherent?(pressure, tg_direction)
          tf_abs[tf] += pressure.abs
          tf_pos[tf] += pressure if pressure > 0
          tf_neg[tf] += pressure.abs if pressure < 0
          tf_contexts[tf] << context
          tf_targets[tf] << tg
        end
      end

      top_coherent = tf_abs.sort_by { |tf, value| [-value, tf] }.first(top_n).collect { |tf, value| "#{tf}(#{value.round(2)})" } * ","
      top_positive = tf_pos.sort_by { |tf, value| [-value, tf] }.first(top_n).collect { |tf, value| "#{tf}(#{value.round(2)})" } * ","
      top_negative = tf_neg.sort_by { |tf, value| [-value, tf] }.first(top_n).collect { |tf, value| "#{tf}(#{value.round(2)})" } * ","
      top_tf_contexts = tf_abs.sort_by { |tf, value| [-value, tf] }.first([top_n, 5].min).collect { |tf, _value| "#{tf}:#{tf_contexts[tf].to_a.sort.first(8) * '|'}" } * ","
      top_tf_targets = tf_abs.sort_by { |tf, value| [-value, tf] }.first([top_n, 5].min).collect { |tf, _value| "#{tf}:#{tf_targets[tf].to_a.sort.first(12) * '|'}" } * ","

      result["M%03d" % mod_i] = [
        entities.length.to_s,
        pos_contexts * ",",
        neg_contexts * ",",
        nonzero_contexts * ",",
        curated_overlap,
        top_coherent,
        top_positive,
        top_negative,
        top_tf_contexts,
        top_tf_targets,
        grn_format_gene_list(entities, 100),
        sig
      ]
    end

    result
  end

end
