module AGS

  dep :change_offsets
  dep :tf_predictions
  task :grns => :array do
    change_offsets = step(:change_offsets).load
    decoupler_jobs = rec_dependencies.select do |dep|
      dep.task_name == :timepoint_decoupler
    end

    self.progress_bar max: decoupler_jobs.length do |bar|
      decoupler_jobs.each do |job|
        treatment, time_point = job.recursive_inputs.values_at :treatment, :time_point

        tfs = {}
        tgs = {}
        edges = {}

        regulome = job.step(:filtered_regulome).path.tsv type: :double
        regulome = regulome.change_key 'source', merge: true
        regulome.unnamed = false

        job.load.each do |tf,activity|
          tfs[tf] = activity.first
          NamedArray.zip_fields(regulome[tf]).each do |tg,weight|
            changes = change_offsets[tg][treatment]
            direction = nil

            # Iterate in reverse and allow matching in previous timepoint if no
            # more recent change is found
            changes.reverse.each do |change|
              dir,_sep, time = change.partition ' '
              time = time.sub 'h', ''
              if time.include? '-'
                start, eend = time.split '-'
                match = true if (start.to_i..eend.to_i).include? time_point
              else
                match = true if time.to_i == time_point
                match = true if time.to_i == time_point - 1
              end

              if match
                direction = dir
                break
              end
            end
            iii [changes, time_point] if direction.nil?
            tgs[tg] = direction
            edges[[tf, tg]*'~'] = weight.to_i
          end
        end


        grn = {tfs: tfs, tgs: tgs, edges: edges}

        target = [treatment, time_point] * '-'
        Open.write file(target + '.json'), grn.to_json
        bar.tick
      end
    end
    files
  end

end
