#!/bin/bash

####Running Gene Flow Simulations Scripts 
set -e

####NOTE: Taxa labeled by their specific epiphet

##Aurulentus
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/aurulentus
#Run Model
gadma -i North-South.sfs -p params_aurulentus.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##Cinerascens
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/cinerascens
#Run Model
gadma -i North-South.sfs -p params_cinerascens.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##Contaminatus
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/contaminatus
#Run Model
gadma -i South-North.sfs -p params_contaminatus.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##desmaresti
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/desmaresti
#Run Model
gadma -i South-North.sfs -p params_desmaresti.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##diadematus
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/diadematus
#Run Model
gadma -i South-North.sfs -p params_diadematus.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##difficilis
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/difficilis
#Run Model
gadma -i South-North.sfs -p params_difficilis.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##lalandi
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/lalandi
#Run Model
gadma -i South-North.sfs -p params_lalandi.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##lateralis
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/lateralis
#Run Model
gadma -i North-South.sfs -p params_lateralis.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##maxillosus
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/maxillosus
#Run Model
gadma -i South-North.sfs -p params_maxillosus.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##obsoletus
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/obsoletus
#Run Model
gadma -i South-North.sfs -p params_obsoletus.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##rubricollis
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/rubricollis
#Run Model
gadma -i South-North.sfs -p params_rubricollis.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##rufosuperciliata
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/rufosuperciliata
#Run Model
gadma -i South-North.sfs -p params_rufosuperciliata.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##thoracica
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/thoracica
#Run Model
gadma -i North-South.sfs -p params_thoracica.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py

##xanthopterus
#Set working directory
cd /work/gjmonaha/miniconda3/envs/gadma_env/3_vcfs/GADMA/xanthopterus
#Run Model
gadma -i North-South.sfs -p params_xanthopterus.txt -o ./results
#Generating Graph
cd ./results
python best_logLL_model_moments_code.py