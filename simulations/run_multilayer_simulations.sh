

observation_depth_L_values=(2)

app_coverage_values=(0.01 0.025 0.005 0.001 0.0005)

fraction_of_private_profiles_values=(0.333333 0.500000 0.66666)

nb_simulations=500


# Compiles the binary file if it does not already exist.
if [[ ! -f bin/multilayer ]]; then
  echo 'compiling bin/multilayer'
  g++ -O3 -std=c++11 multilayer.cpp -o bin/multilayer
fi


# Loops over all available graphs.
networks=(Dartmouth6 Brown11) # UC64 William77 Williams40 Brandeis99 Maine59 UCSC68 Johns_Hopkins55 Vassar85 Vanderbilt48 Duke14 Georgetown15 Rice31 American75 USFCA72 Mich67 Colgate88 Carnegie49 Rochester38 UChicago30 Haverford76 Princeton12 Wesleyan43 Yale4 Caltech36 WashU32 Swarthmore42 Reed98 Simmons81 Bowdoin47 Tulane29 MIT8 Wake73 Pepperdine86 Hamilton46 Bucknell39 Emory27 Vermont70 Trinity100 Santa74 Middlebury45 Wellesley22 Tufts18 Howard90 Oberlin44 Smith60 Amherst41 Villanova62 Lehigh96)
for nname in ${networks[*]}; do
  edgelist_filename=../Facebook100/${nname}.txt.tar.xz
# for edgelist_filename in ../Facebook100/*.txt.tar.xz; do

  network_name=${edgelist_filename##*/}
  network_name=${network_name%%.*}

  output_filename=../results/multilayer/${network_name}.dat

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
  if [[ ! -f ../Facebook100/${network_name}.txt ]]; then
    echo 'extracting edgelist for '${network_name}
    tar xJf ${edgelist_filename}
    mv ${network_name}.txt ../Facebook100/
  fi

  # Runs the script.
  for observation_depth_L in ${observation_depth_L_values[*]}; do
    for app_coverage in ${app_coverage_values[*]}; do
      for fraction_of_private_profiles in ${fraction_of_private_profiles_values[*]}; do
        for adoption_of_passports in 0 0.25 0.50 0.90 0.95 1.00; do

          # When running simulations on a Compute Canada cluster.
          bash launch_multilayer_simulations_on_cedar.sh $network_name $observation_depth_L $app_coverage $fraction_of_private_profiles $adoption_of_passports $nb_simulations $output_filename 0 0.025 1.000000001

          # # When running simulations locally.
          # adoption_within_private_profiles_values=$(seq 0 0.025 1.000000001)
          # for adoption_within_private_profiles in ${adoption_within_private_profiles_values[*]}; do
          #   bin/emergence_of_components ${network_name}.txt $observation_depth_L $app_coverage $fraction_of_private_profiles $adoption_within_private_profiles $adoption_of_passports $nb_simulations $network_name >> $output_filename
          # done

        done
      done
    done
  done

  # # Removes the edgelist filename (running simulations locally).
  # rm ../Facebook100/${network_name}.txt

done
