

observation_depth_L_values=(2)

app_coverage_values=(0.01)

fraction_of_private_profiles_values=(0.333333) # 0.500000 0.800000)

adoption_within_private_profiles_values=$(seq 0 0.025 1.000000001)

nb_simulations=10


# Compiles the binary file if it does not already exist.
if [[ ! -f bin/emergence_of_components ]]; then
  g++ -O3 -std=c++11 src/emergence_of_components.cpp -o bin/emergence_of_components
fi


# output_filename="../../../data_ethics_data/results/emergence_of_components/emergence_of_components.dat"
# # if [[ -f $output_filename ]]; then rm $output_filename; fi


# Loops over all available graphs.
networks=(Dartmouth6 Brown11 UC64 William77 Williams40 Brandeis99 Maine59 UCSC68 Johns_Hopkins55 Vassar85 Vanderbilt48 Duke14 Georgetown15 Rice31 American75 USFCA72 Mich67 Colgate88 Carnegie49 Rochester38 UChicago30 Haverford76 Princeton12 Wesleyan43 Yale4 Caltech36 WashU32 Swarthmore42 Reed98 Simmons81 Bowdoin47 Tulane29 MIT8 Wake73 Pepperdine86 Hamilton46 Bucknell39 Emory27 Vermont70 Trinity100 Santa74 Middlebury45 Wellesley22 Tufts18 Howard90 Oberlin44 Smith60 Amherst41 Villanova62 Lehigh96)
for nname in ${networks[*]}; do
  edgelist_filename=../Facebook100/${nname}.txt.tar.xz
# for edgelist_filename in ../../../data_ethics_data/Facebook100/*.txt.tar.xz; do

  network_name=${edgelist_filename##*/}
  network_name=${network_name%%.*}

  output_filename=results/${network_name}.dat
  # if [[ -f $output_filename ]]; then rm $output_filename; fi

  if [[ ! -f ${network_name}_empty.tmp ]]; then
      touch ${network_name}_empty.tmp
      echo ${network_name}

  # Prints the header of the file if the files does not already exist.
  if [[ ! -f $output_filename ]]; then
    touch $output_filename
    echo -n "#          Name " >> $output_filename
    echo -n "       ObsDepth " >> $output_filename
    echo -n "    AppCoverage " >> $output_filename
    echo -n "   PrivProfFrac " >> $output_filename
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
          ./bin/emergence_of_components ${network_name}.txt $observation_depth_L $app_coverage $fraction_of_private_profiles $adoption_within_private_profiles $nb_simulations $network_name >> $output_filename
        done
      done
    done
  done
  # Removes the edgelist filename
  rm ${network_name}.txt

  fi

done
