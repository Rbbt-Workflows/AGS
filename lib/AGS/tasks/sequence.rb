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
            key = Misc.digest(values)
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
end
