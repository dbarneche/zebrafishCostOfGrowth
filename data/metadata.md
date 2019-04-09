# Datasets for zebrafish cost of growth  
## data/data.csv  
### Variable descriptions  
*tankID:* The ID of the experimental tank, given as temperature_tank number;  
*205A_origin:* The ID of the breeding tank (first month of rearing);  
*Basket:* Basket where fish was kept during the experiment (5 baskets per tankID);  
*fishID:* The combination of columns tankID and Basket;  
*witroxChannel:* Identifier for the Loligo Witrox channel used;  
*respirometryChamber:* Identifier for the acrylic vial used;  
*witroxFileName:* Raw file containing time series respirometry data generated using Loligo Witrox;  
*controlChamber:* ID of control chamber for respirometry assay;  
*Date:* Date of measurement (day/month/year);  
*mass_g:* Measured mass of the individual, in grams;  
*OBS:* Any particular observation.  
  
## data/mitoResp.csv  
### Variable descriptions  
*fishID:* The combination of columns tankID and Basket (matches that of data/data.csv);  
*tempCelsius:* Experimental temperature in Celsius degrees;  
*tankID:* The ID of the experimental tank, given as temperature_tank number (matches that of data/data.csv);  
*beginningState2_O2_umol_l:* Concentration of O~2~ at the beginning of State 2;  
*endState2_O2_umol_l:* Concentration of O~2~ at the end of State 2;  
*state2Time_mins:* Duration of State 2 (in minutes);  
*beginningState3_maxO2_umol_l:* Concentration of maximum O~2~ at the beginning of State 3;  
*endState3_maxO2_umol_l:* Concentration of maximum O~2~ at the end of State 3;  
*maxState3Time_mins:* Duration of maximum State 3 (in minutes);  
*beginningState4_O2_umol_l:* Concentration of O~2~ at the beginning of State 4;  
*endState4_O2_umol_l:* Concentration of O~2~ at the end of State 4;  
*state4Time_mins:* Duration of State 4 (in minutes);  
*startOliO2_umol_l:* Concentration of O~2~ at the beginning of oligomycin addition;  
*endOliO2_umol_l:* Concentration of O~2~ at the end of oligomycin addition;  
*oliTime_mins:* Duration of oligomycin addition (in minutes);  
*O2BeginningS3_O2_umol_l:* Concentration of O~2~ at the beginning of State 3;  
*O2EndS3_O2_umol_l:* Concentration of O~2~ at the beginning of State 3;  
*timeO2DiffusionS3_mins:* Duration of O~2~ diffusion in State 3 (in minutes);  
*proteinSample_mg_ml:* Protein concentration for Bradford assay;  
*protein_per_50ul:* Protein concentration per 50 µL;  
