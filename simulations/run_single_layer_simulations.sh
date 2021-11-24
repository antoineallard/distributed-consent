#!/bin/bash

# Antoine Allard
# antoineallard.info
# Octobre 2021


observation_depth_L_values=(2)

# app_coverage_values=(0.01 0.025 0.005 0.001 0.0005 0.0001)
app_coverage_values=(0.0005 0.0001)

# fraction_of_private_profiles_values=(0.333333 0.500000 0.66666)
fraction_of_private_profiles_values=(0.333333 0.66666)

nb_simulations=1000


# Compiles the binary file if it does not already exist.
if [[ ! -f bin/single_layer ]]; then
  echo 'compiling bin/single_layer'
  g++ -O3 -std=c++11 single_layer.cpp -o bin/single_layer
fi


# Loops over all available graphs.
for edgelist_filename in ../Facebook100/*.txt.tar.xz; do

  network_name=${edgelist_filename##*/}
  network_name=${network_name%%.*}

  output_filename=../results/single_layer/${network_name}.dat

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
  if [[ ! -f ../Facebook100/${network_name}.txt ]]; then
    echo 'extracting edgelist for '${network_name}
    tar xJf ${edgelist_filename}
    mv ${network_name}.txt ../Facebook100/
  fi

  # Runs the script.
  for observation_depth_L in ${observation_depth_L_values[*]}; do
    for app_coverage in ${app_coverage_values[*]}; do
      for fraction_of_private_profiles in ${fraction_of_private_profiles_values[*]}; do

        # # When running simulations locally.
        # adoption_within_private_profiles_values=$(seq 0 0.025 1.000000001)
        # for adoption_within_private_profiles in ${adoption_within_private_profiles_values[*]}; do
        #   bin/single_layer ${network_name}.txt $observation_depth_L $app_coverage $fraction_of_private_profiles $adoption_within_private_profiles $nb_simulations $network_name >> $output_filename
        # done

        # When running simulations on a Compute Canada cluster.
        bash launch_single_layer_simulations_on_cedar.sh $network_name $observation_depth_L $app_coverage $fraction_of_private_profiles $nb_simulations $output_filename 0 0.025 1.000000001

      done
    done
  done

  # # Removes the edgelist filename (running simulations locally).
  # rm ../Facebook100/${network_name}.txt

done
