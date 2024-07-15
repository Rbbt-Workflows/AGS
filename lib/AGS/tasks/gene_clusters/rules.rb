module AGS

  def self.low_rule(fc_current, fc_next, current_min, next_min, next_multiplier)
    fc_current > current_min && fc_next > next_min && fc_current * next_multiplier > fc_next
  end

  def self.mid_rule(fc_current, fc_next, current_min, next_multiplier)
    fc_current > current_min && fc_next > - fc_current * next_multiplier * (fc_current/current_min)**2
  end

  def self.high_rule(fc_current, current_min)
    fc_current > current_min
  end

  def self.change_start(fc_current, fc_next, low_min, mid_min, high_min, next_min, low_multiplier, mid_multiplier)
    low_rule(fc_current, fc_next, low_min, next_min, low_multiplier) ||
      mid_rule(fc_current, fc_next, mid_min, mid_multiplier) ||
      high_rule(fc_current, high_min)
  end

  def self.cluster_rules
    cluster_rules = {}

    cluster_rules["transient increase 1h to 2h"] = <<-EOF.split("\n")
fcs[0] > 0.15 && (fcs_one[1] < - fcs[0]*0.50)
    EOF

    cluster_rules["transient increase 1h to 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("transient increase 1h to 2h")
fcs [0] > 0.15 && fcs_one[1] < 0 && (fcs_two[2] < - fcs[0]*0.80)
fcs [0] > 0.15 && fcs_one[1] > 0 && (fcs_one[2] < - fcs[1]*0.80)
fcs [0] > 0.10 && fcs_one[1] > 0.03 && (fcs_one[2] < - fcs[1]*0.80)
    EOF

    cluster_rules["transient increase 1h to 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("transient increase 1h to 2h")
"break" if clusters.include?("transient increase 1h to 4h")
fcs [0] > 0.15 && fcs_one[1] < 0 && (fcs_three[3] < - fcs[0]*0.80)
fcs [0] > 0.15 && fcs_one[1] > 0 && (fcs_two[3] < - fcs[1]*0.80)
fcs [0] > 0.15 && fcs_one[1] > 0 && fcs_one[2] > 0 && (fcs_one[3] < - fcs[2]*0.80)
fcs [0] > 0.10 && fcs_one[1] > 0.03 && fcs_one[2] < 0  &&(fcs_two[3] < - fcs[1]*0.80) 
fcs [0] > 0.10 && fcs_one[1] > 0.03 && fcs_one[2] > 0 && (fcs_one[3] < - fcs[2]*0.80) 
    EOF

    cluster_rules["transient decrease 1h to 2h"] = <<-EOF.split("\n")
fcs[0] < -0.15 && (fcs_one[1] > - fcs[0]*0.50)
    EOF

    cluster_rules["transient decrease 1h to 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("transient decrease 1h to 2h")
fcs [0] < -0.15 && fcs_one[1] > 0 && (fcs_two[2] > - fcs[0]*0.80)
fcs [0] < -0.15 && fcs_one[1] < 0 && (fcs_one[2] > - fcs[1]*0.80)
fcs [0] < -0.10 && fcs_one[1] < -0.03 && (fcs_one[2] > - fcs[1]*0.80)
    EOF

    cluster_rules["transient decrease 1h to 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("transient decrease 1h to 2h")
"break" if clusters.include?("transient decrease 1h to 4h")
fcs [0] < -0.15 && fcs_one[1] > 0 && (fcs_three[3] < - fcs[0]*0.80)
fcs [0] < -0.15 && fcs_one[1] < 0 && (fcs_two[3] < - fcs[1]*0.80)
fcs [0] < -0.15 && fcs_one[1] < 0 && fcs_one[2] < 0 && (fcs_one[3] < - fcs[2]*0.80)
fcs [0] < -0.10 && fcs_one[1] < -0.03 && fcs_one[2] > 0  &&(fcs_two[3] < - fcs[1]*0.80) 
fcs [0] < -0.10 && fcs_one[1] < -0.03 && fcs_one[2] < 0 && (fcs_one[3] < - fcs[2]*0.80) 
    EOF

    cluster_rules["start increase 1h"] = <<-EOF.split("\n")
fcs[0] > 0.10 && fcs_one[1] > 0.03 && (fcs_one[1] < fcs[0]*8)
fcs[0] > 0.15 && fcs_one[1] > -fcs[0]*0.60*(fcs[0]/0.15)*(fcs[0]/0.15)
fcs [0] > 0.25
    EOF

    cluster_rules["refactored start increase 1h"] = <<-EOF.split("\n")
    AGS.rule(fcs_one[0], fcs_one[1], low_min_1h, mid_min_1h, high_min_1h, next_min_1h, low_multiplier_1h, mid_multiplier_1h)
    EOF

    cluster_rules["start decrease 1h"] = <<-EOF.split("\n")
fcs[0] < -0.10 && fcs_one[1] < -0.03 && (fcs_one[1] > fcs[0]*8)
fcs[0] < -0.15 && fcs_one[1] < -fcs[0]*0.60*(fcs[0]/0.15)*(fcs[0]/0.15)
fcs [0] < -0.25
    EOF

    cluster_rules["refactored start decrease 1h"] = <<-EOF.split("\n")
    AGS.rule(-fcs_one[0], -fcs_one[1], low_min_1h, mid_min_1h, high_min_1h, next_min_1h, low_multiplier_1h, mid_multiplier_1h)
    EOF


    cluster_rules["transient increase 2h to 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("start increase 1h")
fcs_one[1] > 0.20 && (fcs_one[2] < -fcs_one[1]*0.80)
    EOF

    cluster_rules["transient increase 2h to 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("transient increase 2h to 4h")
fcs_one[1] > 0.20 && fcs_one[2] < 0 && (fcs_two[3] < -fcs_one[1]*0.80)
fcs_one[1] > 0.20 && fcs_one[2] > 0 && (fcs_one[3] < -fcs_two[2]*0.80)
fcs_one[1] > 0.10 && fcs_one[2] > 0.10 && (fcs_one[3] < -fcs_two[2]*0.80)
    EOF

    cluster_rules["transient decrease 2h to 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("start decrease 1h")
fcs_one[1] < -0.20 && (fcs_one[2] > -fcs_one[1]*0.80)
    EOF

    cluster_rules["transient decrease 2h to 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("transient decrease 2h to 4h")
fcs_one[1] < -0.20 && fcs_one[2] > 0 && (fcs_two[3] > -fcs_one[1]*0.80)
fcs_one[1] < -0.20 && fcs_one[2] < 0 && (fcs_one[3] > -fcs_two[2]*0.80)
fcs_one[1] < -0.10 && fcs_one[2] < 0.10 && (fcs_one[3] > -fcs_two[2]*0.80)
    EOF

    cluster_rules["start increase 2h"] = <<-EOF.split("\n")
"break" if clusters.include?("start increase 1h")
fcs_one[1] > 0.10 && fcs_one[2] > 0.10 && (fcs_one[2] < fcs_one[1]*8) 
fcs_one[1] > 0.20 && (fcs_one[2] > - fcs_one[1]*0.70*(fcs_one[1]/0.20)*(fcs_one[1]/0.20)) 
fcs_one[1] > 0.30
    EOF

    cluster_rules["refactored start increase 2h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start increase 1h")
    AGS.rule(fcs_one[1], fcs_one[2], low_min_2h, mid_min_2h, high_min_2h, next_min_2h, low_multiplier_2h, mid_multiplier_2h)
    EOF

    cluster_rules["start decrease 2h"] = <<-EOF.split("\n")
"break" if clusters.include?("start decrease 1h")
fcs_one[1] < -0.10 && fcs_one[2] < -0.10 && (fcs_one[2] > fcs_one[3]*8) 
fcs_one[1] < -0.20 && (fcs_one[2] < - fcs_one[1]*0.70*(fcs_one[1]/0.20)*(fcs_one[1]/0.20)) 
fcs_one[1] < -0.30
    EOF

    cluster_rules["refactored start decrease 2h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start decrease 1h")
    AGS.rule(-fcs_one[1], -fcs_one[2], low_min_2h, mid_min_2h, high_min_2h, next_min_2h, low_multiplier_2h, mid_multiplier_2h)
    EOF


    cluster_rules["transient increase 4h to 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start increase 2h")
fcs_one[2] > 0.20 && (fcs_one[3] < -fcs_one[2]*0.80)
    EOF

    cluster_rules["transient decrease 4h to 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start decrease 2h")
fcs_one[2] < -0.20 && (fcs_one[3] > -fcs_one[2]*0.80)
    EOF

    cluster_rules["start increase 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("start increase 1h") && ! clusters.include?("transient increase 1h to 2h") 
"break" if clusters.include?("start increase 2h") && ! clusters.include?("transient increase 1h to 2h")
fcs_one[2] > 0.10 && fcs_one[3] > 0.10 && (fcs_one[4] < fcs_one[3]*10) 
fcs_one[2] > 0.20 && (fcs_one[3] > - fcs_one[2]*0.90*(fcs_one[2]/0.20)*(fcs_one[2]/0.20)) 
fcs_one[2] > 0.30    
    EOF

    cluster_rules["refactored start increase 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start increase 1h") && ! clusters.include?("transient increase 1h to 2h") 
"break" if clusters.include?("refactored start increase 2h") && ! clusters.include?("transient increase 1h to 2h")
    AGS.rule(fcs_one[2], fcs_one[3], low_min_4h, mid_min_4h, high_min_4h, next_min_4h, low_multiplier_4h, mid_multiplier_4h)
    EOF

    cluster_rules["start decrease 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("start decrease 1h") && ! clusters.include?("transient decrease 1h to 2h") 
"break" if clusters.include?("start decrease 2h") && ! clusters.include?("transient decrease 1h to 2h") 
fcs_one[2] < -0.10 && fcs_one[3] < -0.10 && (fcs_one[4] > fcs_one[3]*10) 
fcs_one[2] < -0.20 && (fcs_one[3] < - fcs_one[2]*0.90*(fcs_one[2]/0.20)*(fcs_one[2]/0.20)) 
fcs_one[2] < -0.30
    EOF

    cluster_rules["refactored start decrease 4h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start decrease 1h") && ! clusters.include?("transient decrease 1h to 2h") 
"break" if clusters.include?("refactored start decrease 2h") && ! clusters.include?("transient decrease 1h to 2h")
    AGS.rule(-fcs_one[2], -fcs_one[3], low_min_4h, mid_min_4h, high_min_4h, next_min_4h, low_multiplier_4h, mid_multiplier_4h)
    EOF

    cluster_rules["start increase 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("start increase 1h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 2h to 4h")
"break" if clusters.include?("start increase 2h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 2h to 4h")
"break" if clusters.include?("start increase 4h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 2h to 4h")
fcs_one[3] > 0.10 && fcs_one[4] > 0.15 && (fcs_one[4] < fcs_one[3]*12)
fcs_one[3] > 0.20 && (fcs_one[4] > - fcs_one[3]*0.90*(fcs_one[3]/0.20)*(fcs_one[3]/0.20)) 
fcs_one[3] > 0.30 
    EOF

    cluster_rules["refactored start increase 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start increase 1h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 2h to 4h")
"break" if clusters.include?("refactored start increase 2h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 2h to 4h")
"break" if clusters.include?("refactored start increase 4h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 2h to 4h")
    AGS.rule(fcs_one[3], fcs_one[4], low_min_8h, mid_min_8h, high_min_8h, next_min_8h, low_multiplier_8h, mid_multiplier_8h)
    EOF

    cluster_rules["start decrease 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("start decrease 1h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 2h to 4h")
"break" if clusters.include?("start decrease 2h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 2h to 4h")
"break" if clusters.include?("start decrease 4h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 2h to 4h")
fcs_one[3] < -0.10 && fcs_one[4] < -0.10
fcs_one[3] < -0.20 && (fcs_one[4] < - fcs_one[3]*0.90*(fcs_one[3]/0.20)*(fcs_one[3]/0.20)) 
fcs_one[3] < -0.30
    EOF

    cluster_rules["refactored start decrease 8h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start decrease 1h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 2h to 4h")
"break" if clusters.include?("refactored start decrease 2h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 2h to 4h")
"break" if clusters.include?("refactored start decrease 4h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 2h to 4h")
    AGS.rule(-fcs_one[3], -fcs_one[4], low_min_4h, mid_min_4h, high_min_4h, next_min_4h, low_multiplier_4h, mid_multiplier_4h)
    EOF

    cluster_rules["start increase 24h"] = <<-EOF.split("\n")
"break" if clusters.include?("start increase 1h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
"break" if clusters.include?("start increase 2h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
"break" if clusters.include?("start increase 4h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
"break" if clusters.include?("start increase 8h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
fcs_one[4] > 0.30
    EOF

    cluster_rules["refactored start increase 24h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start increase 1h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
"break" if clusters.include?("refactored start increase 2h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
"break" if clusters.include?("refactored start increase 4h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
"break" if clusters.include?("refactored start increase 8h") && ! clusters.include?("transient increase 1h to 2h") && ! clusters.include?("transient increase 1h to 4h") && ! clusters.include?("transient increase 1h to 8h") && ! clusters.include?("transient increase 2h to 4h") && ! clusters.include?("transient increase 2h to 8h") && ! clusters.include?("transient increase 4h to 8h")
    AGS.rule(fcs_one[4], 0, low_min_24h, mid_min_24h, high_min_24h, next_min_24h, low_multiplier_24h, mid_multiplier_24h)

    EOF


    cluster_rules["start decrease 24h"] = <<-EOF.split("\n")
"break" if clusters.include?("start decrease 1h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")
"break" if clusters.include?("start decrease 2h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")
"break" if clusters.include?("start decrease 4h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")
"break" if clusters.include?("start decrease 8h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")

fcs_one[4] < -0.30
    EOF

    cluster_rules["refactored start decrease 24h"] = <<-EOF.split("\n")
"break" if clusters.include?("refactored start decrease 1h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")
"break" if clusters.include?("refactored start decrease 2h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")
"break" if clusters.include?("refactored start decrease 4h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")
"break" if clusters.include?("refactored start decrease 8h") && ! clusters.include?("transient decrease 1h to 2h") && ! clusters.include?("transient decrease 1h to 4h") && ! clusters.include?("transient decrease 1h to 8h") && ! clusters.include?("transient decrease 2h to 4h") && ! clusters.include?("transient decrease 2h to 8h") && ! clusters.include?("transient decrease 4h to 8h")
    AGS.rule(-fcs_one[4], 0, low_min_24h, mid_min_24h, high_min_24h, next_min_24h, low_multiplier_24h, mid_multiplier_24h)
    EOF

    cluster_rules
  end
end
