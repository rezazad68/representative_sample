% File paths

%% sample from NASCA
folder = 'NASCA';
add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/NASCA/S04_Results_Factor1,4/';
sample = 'S4_TB102+P450_Resazurin_cE-8M_NADPH_cE-7M_1_FL_532nm_250mW_20ms_all.csv';
Percentage_select = 0.7;
Percentage_drop = 0.3;


% Sample from dSTORM
folder = 'dSTORM';
add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/dSTORM/S16_Results_Right_1,4/';
sample = 'S16_1_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_all.csv';
sample = 'S16_1_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_SameArea_all.csv';
sample = 'S16_2_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_all.csv';
sample = 'S16_2_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_SameArea_all.csv';
Percentage_select = 0.80;
Percentage_drop = 0.2;


add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/dSTORM/S17_Results_Right_1,4/';
sample = 'S17_1_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_all.csv';
sample = 'S17_1_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_SameArea_all.csv';
sample = 'S17_2_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_all.csv';
sample = 'S17_2_TB102+P450_dSTORM_FL_640nm_350mW_NoNileRed_OnlyAlexa+Cys_SameArea_all.csv';

add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/dSTORM/S18_Results_Right_1,4/';
sample = 'S18_1_TB102+P450_dSTORM_PAINT_NR_cE-10_640nm_350mW_20ms_RightChannel_all.csv';
sample = 'S18_2_TB102+P450_dSTORM_PAINT_NR_cE-10_640nm_350mW_20ms_RightChannel_all.csv';
sample = 'S18_2_TB102+P450_dSTORM_PAINT_NR_cE-10_640nm_350mW_20ms_RightChannel_SameArea_all.csv';


% PAINT
folder = 'PAINT';
add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/PAINT/PAINT/S13_Results_Factor1,4/';
sample = 'S10_1_TB102+P450_NASCA+NRPAINT_FL_PAINT_532nm_250mW_20ms_all.csv';
sample = 'S12_1_TB102+P450_NRPAINT_FL_PAINT_532nm_250mW_20ms_all.csv';
sample = 'S10_1_TB102+P450_NASCA+NRPAINT_FL_PAINT_532nm_250mW_20ms_all.csv';
sample = 'S12_1_TB102+P450_NRPAINT_FL_PAINT_532nm_250mW_20ms_all.csv';
sample = 'S12_2_TB102+P450_NRPAINT_FL_PAINT_532nm_250mW_20ms_all.csv';
sample = 'S12_3_TB102+P450_NRPAINT_FL_PAINT_532nm_250mW_20ms_all.csv';
sample = 'S13_1_TB102+P450_NRPAINT_FL_PAINT_532nm_250mW_20ms_all.csv';
Percentage_select = 0.50;
Percentage_drop = 0.2;


% % PAINT
% folder = 'NASCA';
% add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/NASCA/S09_2_Results_Factor1,4/';
% sample = 'S4_TB102+P450_Resazurin_cE-8M_NADPH_cE-7M_1_FL_532nm_250mW_20ms_all.csv';
% sample = 'S4_TB102+P450_Resazurin_cE-8M_NADPH_cE-7M_2_FL_532nm_250mW_20ms_all.csv';
% sample = 'S4_TB102+P450_Resazurin_cE-8M_NADPH_cE-7M_3_FL_532nm_250mW_20ms_all.csv';
% sample = 'S6_TB102+P450_Resazurin_cE-8M_NADPH_cE-7M_1_FL_532nm_250mW_20ms_all.csv';
% sample ='S6_TB102+P450_Resazurin_cE-8M_NADPH_cE-7M_2_FL_532nm_250mW_20ms_all.csv';
% sample = 'S9_2_TB102+P450_NASCA+NRPAINT_1_NASCA_532nm_250mW_20ms_all.csv';
% Percentage_select = 0.60;
% Percentage_drop = 0.2;
% 
% 
% % PAINT and dSTORM
% folder = 'PAINT_dSTORM_combo';
% add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/Combined_PAINT+dSTORM/S15_Results_Factor1,4/';
% sample = 'S14_2_TB102+P450_NRPAINT+dSTORM_FL_dSTORM_640_350mW_20ms_RightChannel_Red_all.csv';
% sample = 'S14_1_TB102+P450_NRPAINT+dSTORM_FL_PAINT_532nm_250mW_dSTORM_640_350mW_20ms_all.csv';
% sample = 'S14_1_TB102+P450_NRPAINT+dSTORM_FL_PAINT_532nm_250mW_dSTORM_640_350mW_20ms_all.csv';
% sample = 'S15_1_TB102+P450_NRPAINT+dSTORM_FL_dSTORM_640_350mW_20ms_RightChannel_Red_dSTORM_all.csv';
% sample = 'S15_1_TB102+P450_NRPAINT+dSTORM_FL_PAINT_532_250mW_20ms_LeftChannel_Green_PAINT_all.csv';
% sample = 'S15_2_TB102+P450_NRPAINT+dSTORM_FL_dSTORM_640_350mW_20ms_RightChannel_Red_dSTORM_all.csv';
% sample = 'S15_2_TB102+P450_NRPAINT+dSTORM_FL_PAINT_532_250mW_20ms_LeftChannel_Green_PAINT_all.csv';
% 
% 
% 
% % PAINT and dSTORM
% folder = 'PAINT_NASCA_combo';
% add = 'E:/lfb research/microgel/sample selection/new_version_new_dataset/Combined_NASCA+PAINT/S11_1_Results_Factor1,4/';
% sample = 'S9_TB102+P450_NASCA+NRPAINT_1_NASCA_532nm_250mW_20ms_all.csv';
% sample = 'S9_TB102+P450_NASCA+NRPAINT_1_PAINT_532nm_250mW_20ms_all.csv';
% sample = 'S11_1_TB102+P450_NASCA+NRPAINT_FL__PAINT_532nm_250mW_20ms_2_all.csv';
% sample = 'S11_1_TB102+P450_NASCA+NRPAINT_FL__PAINT_532nm_250mW_20ms_all.csv';
% sample = 'S11_1_TB102+P450_NASCA+NRPAINT_FL_NASCA_532nm_250mW_20ms_all.csv';


Percentage_select = 0.60;
Percentage_drop = 0.2;


% Hyperparameters
use_spherical_information = false;


raw_pc_path = strcat(add, sample);
% raw_pc_path = strcat('sample_data.txt');
save_path_selected = "results/";

[folder_name, ~, ~] = fileparts(add);
parts = split(folder_name, '/');
exp_name = parts{end};

sample = strtok(sample, '.');

% Combine the last part of the add path with the sample name
exp_name = fullfile(folder, exp_name, sample);




