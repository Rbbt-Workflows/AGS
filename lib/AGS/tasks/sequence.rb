module AGS

  desc "Analyze the time-resolved sequence of transcription factor (TF) activations based on TF activities and a regulome mapping TF → [targets, effects]"
  dep :treatment_tfs
  dep :treatment_tf_consistency
  dep :regulome
  task :sequence_old => :array do |treatment_tfs, regulome|
    treatment_tfs = step(:treatment_tfs).load
    treatment_tf_consistency = step(:treatment_tf_consistency).load

    regulome = step(:regulome).path.tsv key_field: "source", fields: %w(target weight), merge: true, type: :double

	# Build a structure: timepoint => TF => activity_zscore
	tf_series = {}
	treatment_tfs.fields.each do |tp|
	  tf_series[tp] = {}
	end
	treatment_tfs.each do |tf, values|
	  values.each_with_index do |val, i|
		tf_series[treatment_tfs.fields[i]][tf] = val.to_f if val and val != ""
	  end
	end

	timepoints = treatment_tfs.fields

	# Output: sequence of events as array of hashes
	tf_events = []

	timepoints.each_with_index do |tp, i|
	  next_tp = timepoints[i+1] if i+1 < timepoints.size
	  tf_series[tp].each do |tf, tf_activity|
		next unless regulome[tf]
		targets, effects = regulome[tf]
		# Only consider strong activation/repression at this time
		next unless tf_activity.abs > 2   # Put a threshold for call
		targets.zip(effects).each do |tgt, effect|
          next if tf === tgt
		  ["same", "next"].each do |when_tp|
			tgt_time = (when_tp == "same" ? tp : next_tp)
			next unless tgt_time && tf_series[tgt_time][tgt]
			tgt_activity = tf_series[tgt_time][tgt]
			# Determine if regulation corresponds with observed direction
            
            if effect.to_f > 0
			  next unless tf_activity > 2 && tgt_activity > 2
            else
			  next unless tf_activity > 2 && tgt_activity < -2
			end

			tf_events << {
			  :timepoint => tp,
			  :source_tf => tf,
			  :source_activity => tf_activity,
			  :target_tf => tgt,
			  :target_timepoint => tgt_time,
			  :target_activity => tgt_activity,
			  :effect => effect,
			  :type => (when_tp == "same" ? "coincident" : "next")
			}
		  end
		end
	  end
	end

	tf_events
  end

  desc "Analyze the time-resolved sequence of TF activations, following stringent temporal causality and self-consistency criteria"
  dep :treatment_tfs
  dep :treatment_tf_consistency, remove_consecutive: false
  dep :regulome
  task :sequence => :tsv do
    treatment = recursive_inputs[:treatment]
    treatment_tfs = step(:treatment_tfs).load
    treatment_tfs.fields = treatment_tfs.fields.collect{|f| f.split('-').last }

    treatment_tf_consistency = step(:treatment_tf_consistency).load

    regulome = step(:regulome).path.tsv(key_field: "source", fields: %w(target weight), merge: true, type: :double)

    tf_series = {}
    treatment_tfs.fields.each { |tp| tf_series[tp] = {} }
    treatment_tfs.each do |tf, values|
      values.each_with_index do |val, i|
        tf_series[treatment_tfs.fields[i]][tf] = val.to_f if val && val != ""
      end
    end

    # Also build a matrix for self-consistent TFs at each timepoint
    sc_tfs = {}
    treatment_tf_consistency.fields.each do |col|
      next unless col =~ /^Consistent at/
      tp = col.sub("Consistent at ","").sub("h","T").sub("T","")
      tp = "T#{tp}"
      sc_tfs[tp] ||= {}
      treatment_tf_consistency.each do |tf, vals|
        idx = treatment_tf_consistency.fields.index(col)
        sc_tfs[tp][tf] = vals[idx].to_i rescue 0
      end
    end

    timepoints = treatment_tfs.fields
    
    # Change activity threshold for sustained etc
    activity_threshold = 2.3

    # Helper to check for 'sustained' activity (strong at Tn AND at least one subsequent time)
    is_sustained = Proc.new do |tf, idx|
      act_now = tf_series[timepoints[idx]][tf]
      act_now = act_now && act_now.abs > activity_threshold
      if idx < timepoints.length - 1
        later = tf_series[timepoints[idx+1]][tf]
        match_sign = later && act_now && (tf_series[timepoints[idx]][tf] > 0) == (tf_series[timepoints[idx+1]][tf] > 0)
        later = later && later.abs > activity_threshold && match_sign
      else
        later = false
      end
      act_now && later
    end

    tf_events = TSV.setup({}, key_field: "ID", fields: ['Source', 'Source timepoint', 'Source activity',
                                                        'Target', 'Target timepoint', 'Target activity',
                                                        'Effect', 'Offset', 'Type', 'Source self-consistent', 'Target self-consistent']
                         )
    timepoints.each_with_index do |tp, i|
      tf_series[tp].each do |tf_a, tf_a_activity|
        next unless regulome[tf_a]
        targets, effects = regulome[tf_a]
        targets.zip(effects).each do |tf_b, effect|
          next if tf_a == tf_b
          effect = effect.to_f

          (0..2).each do |offset|
            j = i + offset
            next if j >= timepoints.size
            tp_b = timepoints[j]
            tf_b_activity = tf_series[tp_b][tf_b]
            next unless tf_b_activity

            # --- Interaction sign plausibility check ---
            if tf_a_activity > 0
              if effect > 0
                next unless tf_a_activity > activity_threshold && tf_b_activity > activity_threshold
              else
                next unless tf_a_activity > activity_threshold && tf_b_activity < -activity_threshold
              end
            else
              if effect > 0
                next unless tf_a_activity < -activity_threshold && tf_b_activity < -activity_threshold
              else
                next unless tf_a_activity < -activity_threshold && tf_b_activity > activity_threshold
              end
            end

            # --- Temporal offset checks (rules on gaps & self-consistency) ---

            # Require TFb to be self-consistent at this offset (Tn,Tn+1,Tn+2)
            sc_ok = (
              (sc_tfs[tp_b] && sc_tfs[tp_b][tf_b] && sc_tfs[tp_b][tf_b] == 1)
            )

            sc_ok_a = (
              (sc_tfs[tp] && sc_tfs[tp][tf_a] && sc_tfs[tp][tf_a] == 1)
            )

            # Prune: Disqualify Tn+2 links if Tn+2==T24
            if offset==2 && (tp_b == "T24" || tp_b == "T24")
              next
            end

            # Prune: Disqualify T24 if SC is only at T1/2/4, and disqualify T8 if SC only at T1/2

            #if tp_b =~ /T24/
            #  early_consistent = (sc_tfs["T1"] && sc_tfs["T1"][tf_b] && sc_tfs["T1"][tf_b] == 1) ||
            #                     (sc_tfs["T2"] && sc_tfs["T2"][tf_b] && sc_tfs["T2"][tf_b] == 1) ||
            #                     (sc_tfs["T4"] && sc_tfs["T4"][tf_b] && sc_tfs["T4"][tf_b] == 1)
            #  next if early_consistent
            #elsif tp_b =~ /T8/
            #  early_consistent = (sc_tfs["T1"] && sc_tfs["T1"][tf_b] && sc_tfs["T1"][tf_b] == 1) ||
            #                     (sc_tfs["T2"] && sc_tfs["T2"][tf_b] && sc_tfs["T2"][tf_b] == 1)
            #  next if early_consistent
            #end

            # If same timepoint (Tn), require 'sustained' activity
            if offset == 0 && !is_sustained.call(tf_b, i)
              next
            end

            type = (offset==0 ? "sustained/coincident" : "delay_#{offset}")

            values = [tf_a, tp, tf_a_activity, tf_b, tp_b, tf_b_activity, effect, offset, type, sc_ok_a, sc_ok]
            key = treatment + '-' + Misc.digest(values)
            tf_events[key] = values
          end
        end
      end
    end
    tf_events
  end
  
  dep :sequence
  dep :change_offsets
  task :sequence_with_changes => :tsv do
    treatment = self.recursive_inputs[:treatment]
    tsv = step(:sequence).load
    offsets = step(:change_offsets).path.tsv fields: [treatment], type: :flat

    tsv.add_field 'Source Changes' do |key,values|
      source = values['Source'].first
      offsets[source]
    end

    tsv.add_field 'Target Changes' do |key,values|
      target = values['Target'].first
      offsets[target]
    end

    tsv
  end

  desc "Extract relevant chained TF activation sequences by filtering pairwise events (eligible pairs, self-consistency, activity cutoff) and walking backward from T24 through T8, T4, T2, T1"
  dep :sequence, jobname: 'Default'
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :exclude_same_timepoint, :boolean, "Exclude T2-T2 and T4-T4 self-sustaining links (consider to not use)", true
  input :include_orphan_targets, :boolean, "Add orphan targets at each stage that have no downstream link", false
  input :include_t1_sources, :boolean, "Include T1 TF sources (optional early layer)", false
  input :broad_mode, :boolean, "Keep all activity values; annotate pass/fail instead of hard cutoff filter", false
  task :chained_sequences => :tsv do |activity_cutoff, exclude_same_timepoint, include_orphan_targets, include_t1_sources, broad_mode|
    require 'set'

    seq = step(:sequence).load

    # --- Eligible temporal pairs ---
    eligible_pairs = Set.new([
      ["T1", "T2"], ["T1", "T4"],
      ["T2", "T2"], ["T2", "T4"],
      ["T4", "T4"], ["T4", "T8"],
      ["T8", "T8"], ["T8", "T24"],
    ])
    if exclude_same_timepoint
      eligible_pairs.subtract([["T2", "T2"], ["T4", "T4"]])
    end

    # Helper to parse self-consistent values (stored as "true"/"false" strings)
    sc_true = lambda { |v| v.to_s == "true" || v.to_s == "1" }

    # --- Phase A: Filter individual pairwise events ---
    filtered = []
    seq.through do |key, values|
      source      = values["Source"].first
      source_tp   = values["Source timepoint"].first
      source_act  = values["Source activity"].first.to_f
      target      = values["Target"].first
      target_tp   = values["Target timepoint"].first
      target_act  = values["Target activity"].first.to_f
      source_sc   = sc_true.call(values["Source self-consistent"].first)
      target_sc   = sc_true.call(values["Target self-consistent"].first)

      # 1. Eligible temporal pair check
      next unless eligible_pairs.include?([source_tp, target_tp])

      # 2. Self-consistency rules
      #    Target must always be self-consistent
      next unless target_sc
      #    Source must be self-consistent for T4 and T8; T1 and T2 may be non-SC
      if source_tp == "T4" || source_tp == "T8"
        next unless source_sc
      end

      # 3. Activity cutoff (hard filter unless broad_mode)
      passes_cutoff = source_act.abs >= activity_cutoff && target_act.abs >= activity_cutoff
      next unless broad_mode || passes_cutoff

      filtered << {
        key: key,
        source: source, source_tp: source_tp, source_act: source_act,
        target: target, target_tp: target_tp, target_act: target_act,
        source_sc: source_sc, target_sc: target_sc,
        effect: values["Effect"].first.to_f,
        offset: values["Offset"].first.to_i,
        type_str: values["Type"].first,
        passes_cutoff: passes_cutoff,
      }
    end

    # --- Phase B: Backward-walk chain assembly ---
    # Walk from T24 backward through T8, T4, T2.
    # At each stage, keep only events whose target TF (at that timepoint)
    # was already discovered as a source in a later stage.
    # T1 sources are handled in an optional final step.

    discovered = Set.new   # Set of [tf, timepoint] pairs in the chain
    kept_keys  = Set.new
    chain_step = {}        # key -> human-readable stage label

    # Separate T1-source events (handled in optional final step)
    non_t1      = filtered.reject { |e| e[:source_tp] == "T1" }
    t1_events   = filtered.select { |e| e[:source_tp] == "T1" }

    # Group non-T1 events by target timepoint
    by_target_tp = {}
    non_t1.each do |e|
      (by_target_tp[e[:target_tp]] ||= []) << e
    end

    # Process stages from latest to earliest target timepoint
    ["T24", "T8", "T4", "T2"].each do |stage_tp|
      events = by_target_tp[stage_tp] || []
      events.each do |e|
        target_pair = [e[:target], stage_tp]
        is_anchor   = (stage_tp == "T24")
        is_connected = is_anchor || discovered.include?(target_pair)

        if is_connected
          kept_keys << e[:key]
          chain_step[e[:key]] = "#{stage_tp}<-#{e[:source_tp]}"
          discovered << [e[:source], e[:source_tp]]
        elsif include_orphan_targets
          # Orphan: target not connected to chain yet, but add anyway
          kept_keys << e[:key]
          chain_step[e[:key]] = "#{stage_tp}<-#{e[:source_tp]} (orphan)"
          discovered << [e[:target], e[:target_tp]]
          discovered << [e[:source], e[:source_tp]]
        end
      end
    end

    # Optional: T1 source layer
    # T1 events have target at T2 or T4 (from eligible pairs).
    # Include only those whose target is already in the discovered chain.
    if include_t1_sources
      t1_events.each do |e|
        target_pair = [e[:target], e[:target_tp]]
        if discovered.include?(target_pair)
          kept_keys << e[:key]
          chain_step[e[:key]] = "#{e[:target_tp]}<-T1"
        end
      end
    end

    # --- Build output TSV ---
    out_fields = seq.fields.dup + ["Chain step"]
    out_fields << "Passes cutoff" if broad_mode

    result = TSV.setup({}, key_field: "ID", fields: out_fields, type: :list)

    filtered.each do |e|
      next unless kept_keys.include?(e[:key])
      row = [
        e[:source], e[:source_tp], e[:source_act],
        e[:target], e[:target_tp], e[:target_act],
        e[:effect], e[:offset], e[:type_str],
        e[:source_sc], e[:target_sc],
        chain_step[e[:key]],
      ]
      row << e[:passes_cutoff] if broad_mode
      result[e[:key]] = row
    end

    result
  end
end
