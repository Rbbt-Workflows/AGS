require 'json'
require 'set'

module AGS

  #{{{ COMMON NETWORK

  desc "Build the union (common) network across all treatments from chained_sequences"
  dep :chained_sequences, treatment: :placeholder do |jobname, options|
    TREATMENTS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  task :common_network => :tsv do
    edges = {}
    dependencies.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double
      tsv.through do |key, values|
        source   = values["Source"].first
        target   = values["Target"].first
        effect   = values["Effect"].first.to_f
        edge_key = [source, target]
        edges[edge_key] ||= { effect: effect, treatments: Set.new }
        edges[edge_key][:treatments] << treatment
      end
    end

    result = TSV.setup({}, key_field: "Edge", fields: ["Source", "Target", "Effect", "Treatments"], type: :list)
    edges.sort.each do |(source, target), info|
      key = "#{source}~#{target}"
      result[key] = [source, target, info[:effect], info[:treatments].to_a.sort * ","]
    end
    result
  end

  #{{{ SIF EXPORT

  desc "Export the common network in SIF (Simple Interaction File) format"
  dep :common_network
  task :network_sif => :text do
    net = step(:common_network).load
    net = net.to_double unless net.type == :double

    lines = []
    net.through do |key, values|
      source = values["Source"].first
      target = values["Target"].first
      effect = values["Effect"].first.to_f
      interaction = effect > 0 ? "activates" : "represses"
      lines << [source, interaction, target] * "\t"
    end

    sif_text = lines.sort * "\n" + "\n"
    Open.write(file("common_network.sif"), sif_text)
    sif_text
  end

  #{{{ LAYOUT

  desc "Compute a fixed node layout for the common network using graphviz (sfdp/neato)"
  input :engine, :select, "Graphviz layout engine", "sfdp", :select_options => %w(sfdp neato fdp circo dot)
  input :scale, :float, "Scale factor for node spacing", 1.0
  dep :network_sif
  task :network_layout => :tsv do |engine, scale|
    net = step(:common_network).load
    net = net.to_double unless net.type == :double

    nodes = Set.new
    edges = []
    net.through do |key, values|
      source = values["Source"].first
      target = values["Target"].first
      effect = values["Effect"].first.to_f
      nodes << source
      nodes << target
      edges << [source, target, effect]
    end

    # Build DOT input
    dot = []
    dot << "digraph G {"
    dot << "  overlap=false; splines=true; rankdir=LR;"
    dot << "  node [shape=ellipse, width=0.8, height=0.4, fixedsize=false];"
    edges.each do |s, t, eff|
      color = eff > 0 ? "#D73027" : "#4575B4"
      dot << %(  "#{s}" -> "#{t}" [color="#{color}", penwidth=0.5];)
    end
    dot << "}"

    dot_input = dot * "\n"
    Open.write(file("layout_input.dot"), dot_input)

    # Run graphviz engine to get plain output with positions
    cmd = "#{engine} -Tplain -Gscale=#{scale}"
    io = CMD.cmd(cmd, :in => dot_input, :pipe => true, :no_fail => true)
    plain_output = io.read
    io.join

    # Parse plain output: node lines look like: node NAME X Y WIDTH HEIGHT ...
    layout = TSV.setup({}, key_field: "Node", fields: ["X", "Y"], type: :list)
    plain_output.each_line do |line|
      parts = line.split
      next unless parts[0] == "node"
      node_name = parts[1]
      x = parts[2].to_f
      y = parts[3].to_f
      layout[node_name] = [x, y]
    end

    Open.write(file("network_layout.tsv"), layout.to_s)
    layout
  end

  #{{{ NETWORK PANELS (per treatment x timepoint SVG)

  desc "Generate one SVG panel per treatment x timepoint showing TF activities on the fixed common network layout"
  input :activity_cutoff, :float, "Minimum absolute activity to show colored node", 0.0
  dep :common_network
  dep :network_layout
  dep :tf_predictions
  dep :chained_sequences, treatment: :placeholder do |jobname, options|
    TREATMENTS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  task :network_panels => :array do |activity_cutoff|
    net     = step(:common_network).load
    net     = net.to_double unless net.type == :double
    layout  = step(:network_layout).load
    preds   = step(:tf_predictions).load

    # Build per-treatment edge sets from chained_sequences deps
    treatment_edges = {}
    chained_deps = dependencies.select { |d| d.task_name == :chained_sequences }
    chained_deps.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      edge_set = Set.new
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double
      tsv.through do |key, values|
        source = values["Source"].first
        target = values["Target"].first
        edge_set << [source, target]
      end
      treatment_edges[treatment] = edge_set
    end

    # Build node list and edges from common network
    all_nodes = Set.new
    all_edges = []
    net.through do |key, values|
      source = values["Source"].first
      target = values["Target"].first
      effect = values["Effect"].first.to_f
      all_nodes << source
      all_nodes << target
      all_edges << [source, target, effect]
    end

    # Collect treatment-timepoint columns from tf_predictions
    columns = preds.fields.select { |f| f =~ /-T\d+$/ }

    svg_files = []

    FIGURE_TREATMENT_ORDER.each do |treatment|
      edge_set = treatment_edges[treatment] || Set.new

      TIME_POINTS.each do |tp|
        col = "#{treatment}-T#{tp}"
        next unless columns.include?(col)

        # Get activities for this treatment-timepoint
        activities = {}
        col_idx = preds.fields.index(col)
        preds.each do |tf, values|
          val = values[col_idx]
          activities[tf] = val.to_f if val && val != ""
        end

        # Build SVG directly with embedded positions
        svg = build_network_svg(layout, all_nodes, all_edges, edge_set, activities, activity_cutoff, treatment, tp)

        filename = "#{treatment}-T#{tp}.svg"
        Open.write(file(filename), svg)
        svg_files << file(filename).to_s
      end
    end

    # Also write a combined index
    Open.write(file("index.txt"), svg_files * "\n" + "\n")

    svg_files
  end

  helper :activity_to_color do |activity, cutoff|
    if activity && activity.abs >= cutoff && activity.abs >= 2.5
      if activity > 0
        intensity = [[activity.abs / 8.0, 1.0].min, 0.3].max
        r = (215 + (127 - 215) * intensity).round
        g = (48 + (0 - 48) * intensity).round
        b = (39 + (0 - 39) * intensity).round
        sprintf("#%02X%02X%02X", r, g, b)
      else
        intensity = [[activity.abs / 8.0, 1.0].min, 0.3].max
        r = (69 + (8 - 69) * intensity).round
        g = (117 + (30 - 117) * intensity).round
        b = (180 + (107 - 180) * intensity).round
        sprintf("#%02X%02X%02X", r, g, b)
      end
    elsif activity && activity.abs >= cutoff
      if activity > 0
        "#FDD0D0"
      else
        "#D0D8F0"
      end
    else
      "#F5F5F5"
    end
  end

  # ScoutCoder: when iterating over a TSV loaded with type :list, accessing values[field_name]
  # returns a String directly (not an Array). For :double type, it returns an Array.
  # Use tsv.to_double or call .first on arrays to be safe across both types.

  helper :build_network_svg do |layout, nodes, edges, active_edges, activities, cutoff, treatment, tp|
    # Determine bounding box from layout (include all nodes for consistent canvas)
    xs = nodes.collect { |n| layout[n] ? layout[n][0].to_f : 0.0 }
    ys = nodes.collect { |n| layout[n] ? layout[n][1].to_f : 0.0 }
    iif xs
    iif ys
    min_x, max_x = xs.minmax
    min_y, max_y = ys.minmax

    # Normalize coordinates to a pixel canvas
    margin = 60
    svg_w  = 500
    svg_h  = 400

    range_x = (max_x - min_x)
    range_x = 1.0 if range_x == 0
    range_y = (max_y - min_y)
    range_y = 1.0 if range_y == 0

    node_r = 12.0
    inner_w = svg_w - 2 * margin
    inner_h = svg_h - 2 * margin

    # Compute screen positions
    pos = {}
    nodes.each do |n|
      lx = layout[n] ? layout[n][0].to_f : 0.0
      ly = layout[n] ? layout[n][1].to_f : 0.0
      sx = margin + ((lx - min_x) / range_x) * inner_w
      sy = margin + ((max_y - ly) / range_y) * inner_h
      pos[n] = [sx, sy]
    end

    # Helper: clip line endpoints to node boundaries
    clip_edge = lambda do |sx, sy, tx, ty, r|
      dx = tx - sx
      dy = ty - sy
      dist = Math.sqrt(dx * dx + dy * dy)
      return [sx, sy, tx, ty] if dist == 0
      ux = dx / dist
      uy = dy / dist
      [sx + ux * r, sy + uy * r, tx - ux * r, ty - uy * r]
    end

    svg = []
    svg << %(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)
    svg << %(<svg xmlns="http://www.w3.org/2000/svg" width="#{svg_w}" height="#{svg_h}" viewBox="0 0 #{svg_w} #{svg_h}">)
    svg << %(<rect width="#{svg_w}" height="#{svg_h}" fill="white"/>)

    # Title
    label = "#{AGS.figure_treatment_label(treatment.to_s)} T#{tp}"
    svg << %(<text x="#{svg_w / 2}" y="20" text-anchor="middle" font-family="Helvetica" font-size="14" font-weight="bold">#{label}</text>)

    # Arrow marker definitions (must come before use)
    svg << %(<defs>)
    svg << %(<marker id="arrow-act" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto"><polygon points="0,0 10,4 0,8" fill="#D73027"/></marker>)
    svg << %(<marker id="arrow-rep" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto"><polygon points="0,0 10,4 0,8" fill="#4575B4"/></marker>)
    svg << %(</defs>)

    # Inactive edges first (draw under active edges)
    edges.each do |source, target, effect|
      next unless pos[source] && pos[target]
      is_active = active_edges.include?([source, target])
      next if is_active
      sx, sy = pos[source]
      tx, ty = pos[target]
      svg << %(<line x1="#{sx.round(1)}" y1="#{sy.round(1)}" x2="#{tx.round(1)}" y2="#{ty.round(1)}" stroke="#E8E8E8" stroke-width="0.5"/>)
    end

    # Active edges
    edges.each do |source, target, effect|
      next unless pos[source] && pos[target]
      is_active = active_edges.include?([source, target])
      next unless is_active
      sx, sy = pos[source]
      tx, ty = pos[target]
      # Clip to node boundaries
      csx, csy, ctx, cty = clip_edge.call(sx, sy, tx, ty, node_r)
      color = effect > 0 ? "#D73027" : "#4575B4"
      marker_id = effect > 0 ? "arrow-act" : "arrow-rep"
      svg << %(<line x1="#{csx.round(1)}" y1="#{csy.round(1)}" x2="#{ctx.round(1)}" y2="#{cty.round(1)}" stroke="#{color}" stroke-width="1.5" marker-end="url(##{marker_id})"/>)
    end

    # Nodes
    nodes.sort.each do |n|
      next unless pos[n]
      sx, sy = pos[n]
      activity = activities[n]
      fill = activity_to_color(activity, cutoff)
      svg << %(<circle cx="#{sx.round(1)}" cy="#{sy.round(1)}" r="#{node_r}" fill="#{fill}" stroke="#333" stroke-width="0.5"/>)
      # Label
      svg << %(<text x="#{sx.round(1)}" y="#{(sy + node_r + 11).round(1)}" text-anchor="middle" font-family="Helvetica" font-size="8">#{n}</text>)
    end

    svg << %(</svg>)
    svg * "\n"
  end

  #{{{ COMBINED SVG GRID

  desc "Combine all per-treatment-timepoint SVG panels into a single grid SVG"
  dep :network_panels
  task :network_grid_svg => :binary do
    panel_files = step(:network_panels).load

    treatments = FIGURE_TREATMENT_ORDER
    timepoints = TIME_POINTS

    panels = {}
    panel_files.each do |path|
      basename = File.basename(path, ".svg")
      if basename =~ /^(.+)-T(\d+)$/
        treatment = $1
        tp = $2.to_i
        panels[[treatment, tp]] = path
      end
    end

    # Panel dimensions (match the ones set in build_network_svg)
    panel_w = 500
    panel_h = 400

    margin = 40
    title_h = 20
    cell_w = panel_w + margin
    cell_h = panel_h + margin + title_h

    total_w = (cell_w * timepoints.length).round
    total_h = (cell_h * treatments.length).round

    grid_parts = []
    grid_parts << %(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)
    grid_parts << %(<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="#{total_w}" height="#{total_h}" viewBox="0 0 #{total_w} #{total_h}">)
    grid_parts << %(<rect width="#{total_w}" height="#{total_h}" fill="white"/>)

    treatments.each_with_index do |treatment, row|
      timepoints.each_with_index do |tp, col|
        x = col * cell_w
        y = row * cell_h
        label = "#{AGS.figure_treatment_label(treatment.to_s)} T#{tp}"

        grid_parts << %(<text x="#{x + cell_w / 2}" y="#{y + 15}" text-anchor="middle" font-family="Helvetica" font-size="14" font-weight="bold">#{label}</text>)

        path = panels[[treatment, tp]]
        if path && File.exist?(path)
          inner_svg = Open.read(path)
          # Extract inner content (everything between <svg ...> and </svg>)
          inner_content = inner_svg
            .sub(/<\?xml[^>]*\?>/, "")
            .sub(/^<svg[^>]*>/, "")
            .sub(/<\/svg>\s*$/, "")

          offset_x = x + margin / 2
          offset_y = y + title_h

          grid_parts << %(<g transform="translate(#{offset_x}, #{offset_y})">)
          grid_parts << inner_content
          grid_parts << %(</g>)
        end
      end
    end

    grid_parts << %(</svg>)

    grid_svg = grid_parts * "\n"
    Open.write(file("network_grid.svg"), grid_svg)
    grid_svg
  end

end
