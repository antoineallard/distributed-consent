#!/bin/bash
# @author: Antoine Allard <antoineallard.info>
# @author: Laurent-Hébert Dufresne

observation_depth_L_values=(2)

#app_coverage_values=(0.01)
#app_coverage_values=(0.0050126)
app_coverage_values=(0.0025)

fraction_of_private_profiles_values=(0.333333) # 0.500000 0.800000)

adoption_within_private_profiles_values=$(seq 0 0.025 1.000000001)
#adoption_within_private_profiles_values=$(seq 0 0.1 1.000000001)

#adoption_of_passports_values=$(seq 0.95 0.05 1.00000001)
#now set below

#nb_simulations=700
nb_simulations=300


# Compiles the binary file if it does not already exist.
if [[ ! -f bin/emergence_of_components ]]; then
  g++ -O3 -std=c++11 src/emergence_of_components.cpp -o bin/emergence_of_components
fi


# output_filename="../../../data_ethics_data/results/emergence_of_components/emergence_of_components.dat"
# # if [[ -f $output_filename ]]; then rm $output_filename; fi


# Loops over all available graphs.
# networks=(Dartmouth6 Brown11 UC64 William77 Williams40 Brandeis99 Maine59 UCSC68 Johns_Hopkins55 Vassar85 Vanderbilt48 Duke14 Georgetown15 Rice31 American75 USFCA72 Mich67 Colgate88 Carnegie49 Rochester38 UChicago30 Haverford76 Princeton12 Wesleyan43 Yale4 Caltech36 WashU32 Swarthmore42 Reed98 Simmons81 Bowdoin47 Tulane29 MIT8 Wake73 Pepperdine86 Hamilton46 Bucknell39 Emory27 Vermont70 Trinity100 Santa74 Middlebury45 Wellesley22 Tufts18 Howard90 Oberlin44 Smith60 Amherst41 Villanova62 Lehigh96)
networks=(American75 Amherst41 Auburn71 BC17 BU10 Baylor93 Berkeley13 Bingham82 Bowdoin47 Brandeis99 Brown11 Bucknell39 Cal65 Caltech36 Carnegie49 Colgate88 Columbia2 Cornell5 Dartmouth6 Duke14 Emory27 FSU53 GWU54 Georgetown15 Hamilton46 Harvard1 Haverford76 Howard90 Indiana69 JMU79 JohnsHopkins55 Lehigh96 MIT8 MSU24 MU78 Maine59 Maryland58 Mich67 Michigan23 Middlebury45 Mississippi66 NYU9 Northeastern19 Northwestern25 NotreDame57 Oberlin44 Oklahoma97 Penn94 Pepperdine86 Princeton12 Reed98 Rice31 Rochester38 Rutgers89 Santa74 Simmons81 Smith60 Stanford3 Swarthmore42 Syracuse56 Temple83 Tennessee95 Texas80 Texas84 Trinity100 Tufts18 Tulane29 UC33 UC61 UC64 UCF52 UCLA26 UCSB37 UCSC68 UCSD34 UChicago30 UConn91 UF21 UGA50 UIllinois20 UMass92 UNC28 UPenn7 USC35 USF51 USFCA72 UVA16 Vanderbilt48 Vassar85 Vermont70 Villanova62 Virginia63 Wake73 WashU32 Wellesley22 Wesleyan43 William77 Williams40 Wisconsin87 Yale4)
# networks=(BU10 FSU53 NYU9 UIllinois20 Northeastern19 Texas84 UC61 Harvard1 Michigan23 Indiana69 Cornell5 Syracuse56 GWU54 UC33 Virginia63 UConn91 UF21 UCSB37 JMU79 Rutgers89 Tennessee95 Mississippi66 Stanford3 Maryland58 Temple83 Columbia2 Texas80 Bingham82 UCSD34 Northwestern25 UMass92 UGA50 Auburn71 UNC28 MSU24 NotreDame57 USC35 Wisconsin87 BC17 Cal65 Oklahoma97 JohnsHopkins55 Berkeley13 USF51 UCLA26 Penn94 MU78 UPenn7 UVA16 UCF52 Baylor93)
for nname in ${networks[*]}; do

  edgelist_filename=../Facebook100/${nname}.txt.tar.xz

  network_name=${edgelist_filename##*/}
  network_name=${network_name%%.*}

  output_filename=results/${network_name}.dat
  # if [[ -f $output_filename ]]; then rm $output_filename; fi

  if [[ ! -f results/${network_name}_empty.tmp ]]; then
      touch results/${network_name}_empty.tmp
      echo ${network_name}

      # Prints the header of the file if the files does not already exist.
      if [[ ! -f $output_filename ]]; then
        touch $output_filename
        echo -n "#          Name " >> $output_filename
        echo -n "       ObsDepth " >> $output_filename
        echo -n "    AppCoverage " >> $output_filename
        echo -n "   PrivProfFrac " >> $output_filename
        echo -n "   PassportAdop " >> $output_filename
        echo -n "   AdoptionRate " >> $output_filename
        echo -n "     NbVertices " >> $output_filename
        echo -n "      ObsNbComp " >> $output_filename
        echo -n "     ObsNbType0 " >> $output_filename
        echo -n "     ObsNbType1 " >> $output_filename
        echo -n "     ObsNbType2 " >> $output_filename
        echo -n "   ObsGCNbType0 " >> $output_filename
        echo -n "   ObsGCNbType1 " >> $output_filename
        echo -n "   ObsGCNbType2 " >> $output_filename
        echo -n "     NObsNbComp " >> $output_filename
        echo -n "    NObsNbType0 " >> $output_filename
        echo -n "    NObsNbType1 " >> $output_filename
        echo -n "    NObsNbType2 " >> $output_filename
        echo -n "  NObsGCNbType0 " >> $output_filename
        echo -n "  NObsGCNbType1 " >> $output_filename
        echo -n "  NObsGCNbType2 " >> $output_filename
        echo "" >> $output_filename
      fi

      # Uncompiles the archive containing the edgelist
      tar xJf ${edgelist_filename}
      # tar xJf ${dirname}/${network_name}.txt.tar.xz
      # Runs the script.
      for observation_depth_L in ${observation_depth_L_values[*]}; do
        for app_coverage in ${app_coverage_values[*]}; do
          for fraction_of_private_profiles in ${fraction_of_private_profiles_values[*]}; do
            for adoption_within_private_profiles in ${adoption_within_private_profiles_values[*]}; do
              for adoption_of_passports in 0.50 0.90 0.95 1.00; do
                ./bin/emergence_of_components ${network_name}.txt $observation_depth_L $app_coverage $fraction_of_private_profiles $adoption_within_private_profiles $adoption_of_passports $nb_simulations $network_name >> $output_filename
              done
            done
          done
        done
      done
      # Removes the edgelist filename
      rm ${network_name}.txt

  fi

done
