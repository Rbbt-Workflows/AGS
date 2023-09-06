module AGS
	def self.cluster_rules
		cluster_rules = {}


		#{{{ START INCREASE RULES

		cluster_rules["start increase 1h"] = <<-EOF.split("\n")

fcs[0] > 0.1 && fcs_one[1] > 0 && fcs_one [1] < 1 && fcs[2] >= 0.25
fcs[0] > 0.15 &&  fcs_one[1] > 0.1
fcs[0] > 0.2 &&  fcs_one[1] > 0
fcs [0] > 0.25 &&  (fcs_one[1] > -fcs[0]*0.2)
fcs [0] > 0.3 &&  (fcs_one[1] > -fcs[0]*0.5)
fcs [0] > 0.35 &&  (fcs_one[1] > -fcs[0]*0.7)
fcs [0] > 0.45 &&  (fcs_one[1] > -fcs[0]*0.9)
fcs [0] > 0.5

		EOF


		cluster_rules["start increase 2h"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")

fcs[0] > 0.04 && fcs_one[1] > 0.08 &&  fcs_one [2] > 0.05 && fcs_one[2] < 1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.05 && fcs_one[3] > 0.05 && fcs_one[2] < 1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.07 && fcs_one[3] > 0 && fcs_one[2] < 1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.1 && fcs_one[2] < 1
fcs[0] > 0.04 && fcs_one[1] > 0.12 && (fcs_one [2] > - fcs_one [1]*0.5)
fcs_one[1] > 0.20 &&  fcs_one [2] > 0
fcs_one[1] > 0.25 &&  (fcs_one [2] > - fcs_one [1]*0.5)
fcs_one[1] > 0.3 &&  (fcs_one [2] > - fcs_one [1]*0.7)
fcs_one[1] > 0.35 &&  (fcs_one [2] > - fcs_one [1]*0.9)
fcs_one[1] > 0.45 &&  (fcs_one [2] > - fcs_one [1])
fcs_one[1] > 0.5

		EOF


		cluster_rules["start increase 4h"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start increase 2h")

fcs_two[2] >= 0.12 && fcs_one[1] > 0.04 && fcs_one[2] > 0.05 &&  fcs_one [3] > 0.05  && fcs_one[3] < 1
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.1 && fcs_one [3] > 0
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.2 &&  (fcs_one [3] > - fcs_one [2]*0.9)
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.25 &&  (fcs_one [3] > - fcs_one [2])
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.35 &&  (fcs_one [3] > - fcs_one [2]*1.3)
fcs_one[1] > 0.04  && fcs_one[2] > 0.15 && fcs_one [3] > 0
fcs_one[1] > 0.04  && fcs_one[2] > 0.25 &&  (fcs_one [3] > - fcs_one [2]*0.9)
fcs_one[1] > 0.04  && fcs_one[2] > 0.35 &&  (fcs_one [3] > - fcs_one [2])
fcs_one[1] > 0.04  && fcs_one[2] > 0.4 &&  (fcs_one [3] > - fcs_one [2]*1.3)
fcs_one[2] > 0.15 &&  fcs_one [3] > 0.1 && fcs_one[3] < 1
fcs_one[2] > 0.2 &&  fcs_one [3] > 0
fcs_one[2] > 0.25 &&  (fcs_one [3] > - fcs_one [2]*0.7)
fcs_one[2] > 0.35 &&  (fcs_one [3] > - fcs_one [2])
fcs_one[2] > 0.45 &&  (fcs_one [3] > - fcs_one [2]*1.3)
fcs_one[2] > 0.5

		EOF



		cluster_rules["start increase 8h"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start increase 4h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" increase to new level 4-24h")

fcs_one[3] > 0.15 && fcs_one [4] > 0.1 && fcs_one [4] < 1
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.1 &&  fcs_one [4] > 0
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.2 &&  (fcs_one [4] > - fcs_one [3]*0.9)
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.25 &&  (fcs_one [4] > - fcs_one [3])
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.35 &&  (fcs_one [4] > - fcs_one [3]*1.3)
fcs_one[2] > 0.04 && fcs_one[3] > 0.15 &&  fcs_one [4] > 0
fcs_one[2] > 0.04 && fcs_one[3] > 0.2 &&  (fcs_one [4] > - fcs_one [3]*0.9)
fcs_one[2] > 0.04 && fcs_one[3] > 0.25 &&  (fcs_one [4] > - fcs_one [3])
fcs_one[2] > 0.04 && fcs_one[3] > 0.4 &&  (fcs_one [4] > - fcs_one [3]*1.3)
fcs_one[3] > 0.20 &&  fcs_one [4] > 0
fcs_one[3] > 0.25 &&  (fcs_one [4] > - fcs_one [3]*0.7)
fcs_one[3] > 0.3 &&  (fcs_one [4] > - fcs_one [3]*0.9)
fcs_one[3] > 0.35 &&  (fcs_one [4] > - fcs_one [3])
fcs_one[3] > 0.4 &&  (fcs_one [4] > - fcs_one [3]*1.3)
fcs_one[3] > 0.5

		EOF




		cluster_rules["start increase 24h"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start increase 4h")
"break" if clusters.include?("start increase 8h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" decrease to new level 8-24h")
"break" if clusters.include?(" increase to new level 4-24h")
"break" if clusters.include?(" increase to new level 8-24h")

fcs_one[3] < -0.05  && fcs_ratio[4] <0 && fcs_one[4] > 0.2
fcs_one[3] > 0.05 && fcs_one[4] > 0.2
fcs_one[4] >= 0.25

        EOF


        #{{{ START DECREASE RULES


        cluster_rules["start decrease 1h"] = <<-EOF.split("\n")

fcs[0] < -0.1 && fcs_one[1] < 0 && fcs_one[1] > -1 &&  fcs[3] <=  -0.25
fcs[0] < -0.15 &&  fcs_one[1] <  -0.1
fcs[0] < -0.2 &&  fcs_one[1] <  0
fcs [0] < -0.25 &&  (fcs_one[1] < -fcs[0]*0.2)
fcs [0] < -0.3 &&  (fcs_one[1] < -fcs[0]*0.5)
fcs [0] < -0.35 &&  (fcs_one[1] < -fcs[0]*0.7)
fcs [0] < -0.45 &&  (fcs_one[1] < -fcs[0]*0.9)
fcs [0] < -0.5

        EOF


        cluster_rules["start decrease 2h"] = <<-EOF.split("\n")

"break" if clusters.include?("start decrease 1h")

fcs[0] < -0.04 && fcs_one[1] < -0.08 &&  fcs_one [2] < -0.05 && fcs_one[3] > -1
fcs_one[1] < -0.12 && fcs_one [2] < -0.05 && fcs_one [3] < -0.05 && fcs_one[3] > -1
fcs_one[1] < -0.12 && fcs_one [2] < -0.07 && fcs_one [3] < 0 && fcs_one[3] > -1
fcs_one[1] < -0.12 && fcs_one [2] < -0.1 && fcs_one[3] > -1
fcs[0] < -0.04 && fcs_one[1] < -0.12 &&  (fcs_one [2] < - fcs_one [1]*0.5)
fcs_one[1] < -0.20 &&  fcs_one [2] < 0
fcs_one[1] < -0.25 &&  (fcs_one [2] < - fcs_one [1]*0.5)
fcs_one[1] < -0.3 &&  (fcs_one [2] < - fcs_one [1]*0.7)
fcs_one[1] < -0.35 &&  (fcs_one [2] < - fcs_one [1]*0.9)
fcs_one[1] < -0.45 &&  (fcs_one [2] < - fcs_one [1])
fcs_one[1] < -0.5

        EOF

        cluster_rules["start decrease 4h"] = <<-EOF.split("\n")

"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start decrease 2h")

fcs_two[2] <= -0.12 && fcs_one[1] < -0.04 && fcs_one[2] < -0.05 &&  fcs_one [3] < -0.05 && fcs_one[3] > -1
fcs_one[1] > 0.04  && fcs_one[2] < -0.1 &&  fcs_one [3] <  0
fcs_one[1] > 0.04  && fcs_one[2] < -0.2 &&  (fcs_one [3] < - fcs_one [2]*0.9)
fcs_one[1] > 0.04  && fcs_one[2] < -0.25 &&  (fcs_one [3] < - fcs_one [2])
fcs_one[1] > 0.04  && fcs_one[2] < -0.35 &&  (fcs_one [3] < - fcs_one [2]*1.3)
fcs_one[1] < -0.04 && fcs_one[2] < -0.15 &&  fcs_one [3] <  0
fcs_one[1] < -0.04 && fcs_one[2] < -0.25 &&  (fcs_one [3] < - fcs_one [2]*0.9)
fcs_one[1] < -0.04 && fcs_one[2] < -0.35 &&  (fcs_one [3] < - fcs_one [2])
fcs_one[1] < -0.04 && fcs_one[2] < -0.4 &&  (fcs_one [3] < - fcs_one [2]*1.3)
fcs_one[2] < -0.15 && fcs_one [3] < -0.1 && fcs_one[3] > -1
fcs_one[2] < -0.20 &&  fcs_one [3] <  0
fcs_one[2] < -0.25 &&  (fcs_one [3] < - fcs_one [2]*0.7)
fcs_one[2] < -0.35 &&  (fcs_one [3] < - fcs_one [2])
fcs_one[2] < -0.45 &&  (fcs_one [3] < - fcs_one [2]*1.3)
fcs_one[2] < -0.5

        EOF



        cluster_rules["start decrease 8h"] = <<-EOF.split("\n")

"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start decrease 2h")
"break" if clusters.include?("start decrease 4h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" increase to new level 4-24h")

fcs_one[3] < -0.15 && fcs_one[3] < -0.1 && fcs_one[4] > -1
fcs_one[2] > 0.04 && fcs_one[3] < -0.1 &&  fcs_one[4] <  0
fcs_one[2] > 0.04 && fcs_one[3] < -0.2 &&  (fcs_one[4] < - fcs_one[3]*0.9)
fcs_one[2] > 0.04 && fcs_one[3] < -0.25 &&  (fcs_one[4] < - fcs_one[3])
fcs_one[2] > 0.04 && fcs_one[3] < -0.4 && (fcs_one[4] < - fcs_one[3]*1.3)
fcs_one[2] < -0.04 && fcs_one[3] < -0.15 &&  fcs_one[4] <  0
fcs_one[2] < -0.04 && fcs_one[3] < -0.2 &&  (fcs_one[4] < - fcs_one[3]*0.9)
fcs_one[2] < -0.04 && fcs_one[3] < -0.25 &&  (fcs_one[4] < - fcs_one[3])
fcs_one[2] < -0.04 && fcs_one[3] < -0.4 && (fcs_one[4] < - fcs_one[3]*1.3)
fcs_one[3] < -0.20 &&  fcs_one[4] <  0
fcs_one[3] < -0.25 &&  (fcs_one[4] < - fcs_one[3]*0.7)
fcs_one[3] < -0.3 &&  (fcs_one[4] < - fcs_one[3]*0.9)
fcs_one[3] < -0.35 &&  (fcs_one[4] < - fcs_one[3])
fcs_one[3] < -0.45 &&  (fcs_one[4] < - fcs_one[3]*1.3)
fcs_two[3] <= -0.5

        EOF




        cluster_rules["start decrease 24h"] = <<-EOF.split("\n")

"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start decrease 2h")
"break" if clusters.include?("start decrease 4h")
"break" if clusters.include?("start decrease 8h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" decrease to new level 8-24h")
"break" if clusters.include?(" increase to new level 4-24h")
"break" if clusters.include?(" increase to new level 8-24h")

fcs_one[3] > 0.05 && fcs_ratio[4] < 0 && fcs_one[4] < -0.2
fcs_one[3] < -0.05 && fcs_one[4] < -0.2
fcs_one[4] < -0.25

        EOF
        
				#{{{ ONSET UP RULES

        cluster_rules["onset 2h up"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")

fcs[0] > 0.04 && fcs_one[1] > 0.08 &&  fcs_one [2] > 0.05 && fcs_one[2] < 1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.05 && fcs_one[3] > 0.05 && fcs_one[2] < 1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.07 && fcs_one[3] > 0 && fcs_one[2] < 1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.1 && fcs_one[2] < 1
fcs[0] > 0.04 && fcs_one[1] > 0.12 && (fcs_one [2] > - fcs_one [1]*0.5)
fcs_one[1] > 0.20 &&  fcs_one [2] > 0
fcs_one[1] > 0.25 &&  (fcs_one [2] > - fcs_one [1]*0.5)
fcs_one[1] > 0.3 &&  (fcs_one [2] > - fcs_one [1]*0.7)
fcs_one[1] > 0.35 &&  (fcs_one [2] > - fcs_one [1]*0.9)
fcs_one[1] > 0.45 &&  (fcs_one [2] > - fcs_one [1])
fcs_one[1] > 0.5

        EOF

        cluster_rules["onset 4h up"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start decrease 2h")

fcs_two[2] >= 0.12 && fcs_one[1] > 0.04 && fcs_one[2] > 0.05 &&  fcs_one [3] > 0.05  && fcs_one[3] < 1
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.1 && fcs_one [3] > 0
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.2 &&  (fcs_one [3] > - fcs_one [2]*0.9)
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.25 &&  (fcs_one [3] > - fcs_one [2])
fcs_one[1] < -0.04 && fcs_ratio[2] < 0 && fcs_one[2] > 0.35 &&  (fcs_one [3] > - fcs_one [2]*1.3)
fcs_one[1] > 0.04  && fcs_one[2] > 0.15 && fcs_one [3] > 0
fcs_one[1] > 0.04  && fcs_one[2] > 0.25 &&  (fcs_one [3] > - fcs_one [2]*0.9)
fcs_one[1] > 0.04  && fcs_one[2] > 0.35 &&  (fcs_one [3] > - fcs_one [2])
fcs_one[1] > 0.04  && fcs_one[2] > 0.4 &&  (fcs_one [3] > - fcs_one [2]*1.3)
fcs_one[2] > 0.15 &&  fcs_one [3] > 0.1 && fcs_one[3] < 1
fcs_one[2] > 0.2 &&  fcs_one [3] > 0
fcs_one[2] > 0.25 &&  (fcs_one [3] > - fcs_one [2]*0.7)
fcs_one[2] > 0.35 &&  (fcs_one [3] > - fcs_one [2])
fcs_one[2] > 0.45 &&  (fcs_one [3] > - fcs_one [2]*1.3)
fcs_one[2] > 0.5

        EOF

        cluster_rules["onset 8h up"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start decrease 2h")
"break" if clusters.include?("start increase 4h")
"break" if clusters.include?("start decrease 4h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" increase to new level 4-24h")

fcs_one[3] > 0.15 && fcs_one [4] > 0.1 && fcs_one [4] < 1
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.1 &&  fcs_one [4] > 0
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.2 &&  (fcs_one [4] > - fcs_one [3]*0.9)
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.25 &&  (fcs_one [4] > - fcs_one [3])
fcs_one[2] < -0.04 && fcs_ratio[3] < 0 && fcs_one[3] > 0.35 &&  (fcs_one [4] > - fcs_one [3]*1.3)
fcs_one[2] > 0.04 && fcs_one[3] > 0.15 &&  fcs_one [4] > 0
fcs_one[2] > 0.04 && fcs_one[3] > 0.2 &&  (fcs_one [4] > - fcs_one [3]*0.9)
fcs_one[2] > 0.04 && fcs_one[3] > 0.25 &&  (fcs_one [4] > - fcs_one [3])
fcs_one[2] > 0.04 && fcs_one[3] > 0.4 &&  (fcs_one [4] > - fcs_one [3]*1.3)
fcs_one[3] > 0.20 &&  fcs_one [4] > 0
fcs_one[3] > 0.25 &&  (fcs_one [4] > - fcs_one [3]*0.7)
fcs_one[3] > 0.3 &&  (fcs_one [4] > - fcs_one [3]*0.9)
fcs_one[3] > 0.35 &&  (fcs_one [4] > - fcs_one [3])
fcs_one[3] > 0.4 &&  (fcs_one [4] > - fcs_one [3]*1.3)
fcs_one[3] > 0.5

        EOF

        cluster_rules["onset 24h up"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start decrease 2h")
"break" if clusters.include?("start increase 4h")
"break" if clusters.include?("start decrease 4h")
"break" if clusters.include?("start increase 8h")
"break" if clusters.include?("start decrease 8h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" decrease to new level 8-24h")
"break" if clusters.include?(" increase to new level 4-24h")
"break" if clusters.include?(" increase to new level 8-24h")

fcs_one[3] < -0.05  && fcs_ratio[4] <0 && fcs_one[4] > 0.2
fcs_one[3] > 0.05 && fcs_one[4] > 0.2
fcs_one[4] >= 0.25

        EOF

        #{{{ ONSET DOWN RULES

        cluster_rules["onset 2h down"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")

fcs[0] < -0.04 && fcs_one[1] < -0.08 &&  fcs_one [2] < -0.05 && fcs_one[3] > -1
fcs_one[1] < -0.12 && fcs_one [2] < -0.05 && fcs_one [3] < -0.05 && fcs_one[3] > -1
fcs_one[1] < -0.12 && fcs_one [2] < -0.07 && fcs_one [3] < 0 && fcs_one[3] > -1
fcs_one[1] < -0.12 && fcs_one [2] < -0.1 && fcs_one[3] > -1
fcs[0] < -0.04 && fcs_one[1] < -0.12 &&  (fcs_one [2] < - fcs_one [1]*0.5)
fcs_one[1] < -0.20 &&  fcs_one [2] < 0
fcs_one[1] < -0.25 &&  (fcs_one [2] < - fcs_one [1]*0.5)
fcs_one[1] < -0.3 &&  (fcs_one [2] < - fcs_one [1]*0.7)
fcs_one[1] < -0.35 &&  (fcs_one [2] < - fcs_one [1]*0.9)
fcs_one[1] < -0.45 &&  (fcs_one [2] < - fcs_one [1])
fcs_one[1] < -0.5

        EOF

        cluster_rules["onset 4h down"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start decrease 2h")

fcs_two[2] <= -0.12 && fcs_one[1] < -0.04 && fcs_one[2] < -0.05 &&  fcs_one [3] < -0.05 && fcs_one[3] > -1
fcs_one[1] > 0.04  && fcs_one[2] < -0.1 &&  fcs_one [3] <  0
fcs_one[1] > 0.04  && fcs_one[2] < -0.2 &&  (fcs_one [3] < - fcs_one [2]*0.9)
fcs_one[1] > 0.04  && fcs_one[2] < -0.25 &&  (fcs_one [3] < - fcs_one [2])
fcs_one[1] > 0.04  && fcs_one[2] < -0.35 &&  (fcs_one [3] < - fcs_one [2]*1.3)
fcs_one[1] < -0.04 && fcs_one[2] < -0.15 &&  fcs_one [3] <  0
fcs_one[1] < -0.04 && fcs_one[2] < -0.25 &&  (fcs_one [3] < - fcs_one [2]*0.9)
fcs_one[1] < -0.04 && fcs_one[2] < -0.35 &&  (fcs_one [3] < - fcs_one [2])
fcs_one[1] < -0.04 && fcs_one[2] < -0.4 &&  (fcs_one [3] < - fcs_one [2]*1.3)
fcs_one[2] < -0.15 && fcs_one [3] < -0.1 && fcs_one[3] > -1
fcs_one[2] < -0.20 &&  fcs_one [3] <  0
fcs_one[2] < -0.25 &&  (fcs_one [3] < - fcs_one [2]*0.7)
fcs_one[2] < -0.35 &&  (fcs_one [3] < - fcs_one [2])
fcs_one[2] < -0.45 &&  (fcs_one [3] < - fcs_one [2]*1.3)
fcs_one[2] < -0.5

        EOF

        cluster_rules["onset 8h down"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start decrease 2h")
"break" if clusters.include?("start increase 4h")
"break" if clusters.include?("start decrease 4h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" increase to new level 4-24h")

fcs_one[3] < -0.15 && fcs_one[3] < -0.1 && fcs_one[4] > -1
fcs_one[2] > 0.04 && fcs_one[3] < -0.1 &&  fcs_one[4] <  0
fcs_one[2] > 0.04 && fcs_one[3] < -0.2 &&  (fcs_one[4] < - fcs_one[3]*0.9)
fcs_one[2] > 0.04 && fcs_one[3] < -0.25 &&  (fcs_one[4] < - fcs_one[3])
fcs_one[2] > 0.04 && fcs_one[3] < -0.4 && (fcs_one[4] < - fcs_one[3]*1.3)
fcs_one[2] < -0.04 && fcs_one[3] < -0.15 &&  fcs_one[4] <  0
fcs_one[2] < -0.04 && fcs_one[3] < -0.2 &&  (fcs_one[4] < - fcs_one[3]*0.9)
fcs_one[2] < -0.04 && fcs_one[3] < -0.25 &&  (fcs_one[4] < - fcs_one[3])
fcs_one[2] < -0.04 && fcs_one[3] < -0.4 && (fcs_one[4] < - fcs_one[3]*1.3)
fcs_one[3] < -0.20 &&  fcs_one[4] <  0
fcs_one[3] < -0.25 &&  (fcs_one[4] < - fcs_one[3]*0.7)
fcs_one[3] < -0.3 &&  (fcs_one[4] < - fcs_one[3]*0.9)
fcs_one[3] < -0.35 &&  (fcs_one[4] < - fcs_one[3])
fcs_one[3] < -0.45 &&  (fcs_one[4] < - fcs_one[3]*1.3)
fcs_two[3] <= -0.5

        EOF

        cluster_rules["onset 24h down"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?("start increase 2h")
"break" if clusters.include?("start decrease 2h")
"break" if clusters.include?("start increase 4h")
"break" if clusters.include?("start decrease 4h")
"break" if clusters.include?("start increase 8h")
"break" if clusters.include?("start decrease 8h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" decrease to new level 8-24h")
"break" if clusters.include?(" increase to new level 4-24h")
"break" if clusters.include?(" increase to new level 8-24h")

fcs_one[3] > 0.05 && fcs_ratio[4] < 0 && fcs_one[4] < -0.2
fcs_one[3] < -0.05 && fcs_one[4] < -0.2
fcs_one[4] < -0.25

        EOF

        #{{{ ONSET-OFFSET RULES

        cluster_rules["onset 2h up offset 4h"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?(" decrease to new level 4-24h")
"break" if clusters.include?(" increase to new level 4-24h")

fcs[0] > 0.04 && fcs_one[1] > 0.08 &&  fcs_one [2] > 0.05 && fcs_one[2] < 1 && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.05 && fcs_one[3] > 0.05 && fcs_one[2] < 1 && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.07 && fcs_one[3] > 0 && fcs_one[2] < 1 && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.1 && fcs_one[2] < 1 && fcs_ratio[2] < 0 && fcs_one[3] < -0.1
fcs[0] > 0.04 && fcs_one[1] > 0.12 && (fcs_one [2] > - fcs_one [1]*0.5)  && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.20 &&  fcs_one [2] > 0 && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.25 &&  (fcs_one [2] > - fcs_one [1]*0.5) && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.3 &&  (fcs_one [2] > - fcs_one [1]*0.7) && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.35 &&  (fcs_one [2] > - fcs_one [1]*0.9) && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.45 &&  (fcs_one [2] > - fcs_one [1]) && fcs_ratio[2] < 0 && fcs_one[3] < -0.05
fcs_one[1] > 0.5 && fcs_ratio[2] < 0 && fcs_one[3] < -0.05

        EOF

        cluster_rules["onset 2h up offset 8h"] = <<-EOF.split("\n")

"break" if clusters.include?("start increase 1h")
"break" if clusters.include?("start decrease 1h")
"break" if clusters.include?(" decrease to new level 8-24h")
"break" if clusters.include?(" increase to new level 8-24h")

fcs[0] > 0.04 && fcs_one[1] > 0.08 &&  fcs_one [2] > 0.05 && fcs_one[2] < 1 && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.05 && fcs_one[3] > 0.05 && fcs_one[2] < 1 && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.07 && fcs_one[3] > 0 && fcs_one[2] < 1 && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.12 &&  fcs_one [2] > 0.1 && fcs_one[2] < 1 && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs[0] > 0.04 && fcs_one[1] > 0.12 && (fcs_one [2] > - fcs_one [1]*0.5)   && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.20 &&  fcs_one [2] > 0 && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.25 &&  (fcs_one [2] > - fcs_one [1]*0.5) && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.3 &&  (fcs_one [2] > - fcs_one [1]*0.7) && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.35 &&  (fcs_one [2] > - fcs_one [1]*0.9) && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.45 &&  (fcs_one [2] > - fcs_one [1]) && fcs_ratio[3] < 0 && fcs_one[3] < -0.1
fcs_one[1] > 0.5 && fcs_ratio[3] < 0 && fcs_one[3] < -0.1

        EOF
        
				#{{{ OTHER RULES

        cluster_rules["increase to new level 4-24h"] = <<-EOF.split("\n")

fcs_one[2] > 0.1 && fcs_one[2] < 0.8 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[2] > 0.1 && fcs_one[2] < 0.8 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)
fcs_one[2] > 0.1 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs_ratio[4] > -1.5 && fcs_ratio[4] < 0.5 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[2] > 0.1 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs_ratio[4] > -1.5 && fcs_ratio[4] < 0.5 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)
fcs_one[1] > 0.1 && fcs_one[2] > 0.05 && fcs_one[2] < 0.8 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[1] > 0.1 && fcs_one[2] > 0.05 && fcs_one[2] < 0.8 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)
fcs_one[1] > 0.1 && fcs_one[2] > 0.05 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs_ratio[4] > -1.5 && fcs_ratio[4] < 0.5 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[1] > 0.1 && fcs_one[2] > 0.05 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs_ratio[4] > -1.5 && fcs_ratio[4] < 0.5 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)

        EOF




        cluster_rules["decrease to new level 4-24h"] = <<-EOF.split("\n")

fcs_one[2] < -0.1 && fcs_one[2] > -0.8 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[2] < -0.1 && fcs_one[2] > -0.8 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)
fcs_one[2] < -0.1 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[2] < -0.1 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)
fcs_one[1] < -0.1 && fcs_one[2] < 0 && fcs_one[2] > -0.8 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[1] < -0.1 && fcs_one[2] < 0 && fcs_one[2] > -0.8 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.6 && fcs_ratio[3] < 0.2 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)
fcs_one[1] < -0.1 && fcs_one[2] < 0 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs[3] > 0 && (fcs[4] > fcs[3]*0.7) && (fcs[4] < fcs[3]*1.3)
fcs_one[1] < -0.1 && fcs_one[2] < 0 && fcs_ratio[2] >= 0.2 && fcs_ratio[3] > -0.4 && fcs_ratio[3] < 0.2 && fcs[3] < 0 && (fcs[4] < fcs[3]*0.7) && (fcs[4] > fcs[3]*1.3)

        EOF







        cluster_rules["increase to new level 8-24h"] = <<-EOF.split("\n")

"break" if clusters.include?("increase to new level 4-24h")

fcs_one[3] > 0.1 && fcs_ratio[3] > 0 && fcs_ratio[4] > -0.6 && fcs_ratio[4] < 0.2
fcs_one[2] > 0.1 && fcs_one[3] > 0.05 && fcs_ratio[3] > 0 && fcs_ratio[4] > -0.6 && fcs_ratio[4] < 0.2

        EOF



        cluster_rules["decrease to new level 8-24h"] = <<-EOF.split("\n")

"break" if clusters.include?("decrease to new level 4-24h")

fcs_one[3] < -0.1 && fcs_ratio[3] > 0 && fcs_ratio[4] > -0.6 && fcs_ratio[4] < 0.2
fcs_one[2] < -0.1 && fcs_one[3] < -0.05 && fcs_ratio[3] > 0 && fcs_ratio[4] > -0.6 && fcs_ratio[4] < 0.2

        EOF


        cluster_rules["persistent increase 1-24h"] = <<-EOF.split("\n")

fcs[0] > 0.05 && fcs_one[1] > 0.05 &&  fcs_one [2] > 0.05 && fcs_one[3] > 0 && fcs_one[4] > 0.1 && fcs_one[3] < 1
fcs[0] > 0.05 && fcs_one[1] > 0.05 &&  fcs_one [2] > 0 && fcs_one[3] > 0.05 && fcs_one[4] > 0.1 && fcs_one[3] < 1
fcs_one[1] > 0.05 && fcs_one[2] > 0.05 && fcs_one[3] > 0.05 && fcs_one[4] > 0.1 && fcs_one[3] < 1
fcs_one[1] >= 0.1 && fcs_one[2] > 0.05 && fcs_one[3] > 0 && fcs_one[4] > 0.1
fcs_one[1] >= 0.15 && fcs_one[2] > 0 && fcs_one[3] > 0 && fcs_one[4] > 0.1

        EOF


        cluster_rules["persistent decrease 1-24h"] = <<-EOF.split("\n")

fcs[0] < -0.05 && fcs_one[1] < -0.05 &&  fcs_one [2] < -0.05 && fcs_one[3] < 0 && fcs_ratio[4] > 0.1 && fcs_one[3] > -1
fcs[0] < -0.05 && fcs_one[1] < -0.05 &&  fcs_one [2] < 0 && fcs_one[3] < -0.05 && fcs_ratio[4] > 0.1 && fcs_one[3] > -1
fcs_one[1] <= -0.05 && fcs_one[2] < -0.05 && fcs_one[3] <0 && fcs_ratio[4] > 0.1 && fcs_one[3] > -1
fcs_one[1] <= -0.1 && fcs_one[2] < -0.05 && fcs_one[3] < 0 && fcs_ratio[4] > 0.1
fcs_one[1] <= -0.15 && fcs_one[2] < 0 && fcs_one[3] < 0 && fcs_ratio[4] > 0.1

        EOF



        cluster_rules["persistent increase 2-24h"] = <<-EOF.split("\n")

"break" if clusters.include?("persistent increase 1-24h")

fcs_two[2] >= 0.12 && fcs_one[1] > 0.05 && fcs_one[2] > 0.05 &&  fcs_one [3] > 0  && fcs_ratio[4] > 0.1 && fcs_one[3] < 0.8
fcs_two[2] >= 0.12 && fcs_one[1] > 0.05 && fcs_one[2] > 0 &&  fcs_one [3] > 0.05  && fcs_ratio[4] > 0.1 && fcs_one[3] < 0.8
fcs_one[2] > 0.1 && fcs_one[2] < 0.2 && fcs_one [1] < 0 && fcs_ratio[2] < -0.2 && fcs_one [3] > 0 && fcs_ratio[4] > 0.1 && fcs_one[3] < 0.6
fcs_one[2] >= 0.1 && fcs_one[3] > 0.05 && fcs_ratio[4] > 0.1
fcs_one[2] >= 0.15 && fcs_one[3] > 0 && fcs_ratio[4] > 0.1

        EOF


        cluster_rules["persistent decrease 2-24h"] = <<-EOF.split("\n")

"break" if clusters.include?("persistent decrease 1-24h")

fcs_two[2] <= -0.12 && fcs_one[1] < -0.05 && fcs_one[2] < -0.05 &&  fcs_one [3] < 0 && fcs_ratio[4] > 0.1 && fcs_one[3] > -0.8
fcs_two[2] <= -0.12 && fcs_one[1] < -0.05 && fcs_one[2] < 0 &&  fcs_one [3] < -0.05 && fcs_ratio[4] > 0.1 && fcs_one[3] > -0.8
fcs_one[2] < -0.1 && fcs_one[2] > -0.2 && fcs_one [1] > 0 && fcs_ratio[2] < -0.2 && fcs_one [3] < 0 && fcs_one[3] > -0.8
fcs_one[2] <= -0.1 && fcs_one[3] < -0.05 && fcs_ratio[4] > 0.1
fcs_one[2] <= -0.15 && fcs_one[3] < 0 && fcs_ratio[4] > 0.1

        EOF

        cluster_rules
    end

end
