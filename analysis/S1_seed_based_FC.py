# -*- coding: utf-8 -*-
# @Author  : zhangwei
# @Time    : 2023/6/25 PM 4:01
# @Func    : S1_seed_based_FC

from utils.utils import *

if __name__ == '__main__':
	# Define parameters to be modified
	SUB_ROI_Flag = 'Atlas_xx'  # Replace with the specific subregion name
	CLINICAL_INFO_PATH = '/path/to/clinical_info.csv'  # Replace with the path to the clinical information CSV file
	SAVE_DIR = '/path/to/save/FC_Analysis'  # Replace with the directory path to save the FC analysis results
	BASE_DATA_DIR = '/path/to/base/data'  # Replace with the base directory path of the data
	SEED_MASK_DIR = '/path/to/seed/masks'  # Replace with the directory path of the seed mask files
	FS6_N_VERTEX = 40962  # Number of vertices for fsaverage6
	LH_FS6_N_VERTEX = 13654 * 3  # Number of vertices for lh_fs6
	RH_FS6_N_VERTEX = 13654 * 3  # Number of vertices for rh_fs6
	SEED_LIST = ['GPe', 'GPi', 'Pu', 'STN', 'VIM', 'SN']  # List of seeds

	# Load clinical information
	clinical_info = pd.read_csv(CLINICAL_INFO_PATH)
	subject_list = clinical_info['Subject_ID'].tolist()

	# e.g. DBSControl dataset
	flags = ['DBSControl']
	for flag in flags:
		flag_save_path = os.path.join(SAVE_DIR, f'FC_{flag}_{SUB_ROI_Flag}')
		if not os.path.exists(flag_save_path):
			os.makedirs(flag_save_path)

		for indx, sub in enumerate(subject_list):
			print(sub)
			# Load surface BOLD data
			rh_surf_bold_files = sorted(glob(os.path.join(BASE_DATA_DIR, sub, '*', 'Preprocess', 'surf',
			                                              'rh.*_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6.nii.gz')))
			rh_surf_bold = get_surf_concate(rh_surf_bold_files, RH_FS6_N_VERTEX)

			lh_surf_bold_files = sorted(glob(os.path.join(BASE_DATA_DIR, sub, '*', 'Preprocess', 'surf',
			                                              'lh.*_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6.nii.gz')))
			lh_surf_bold = get_surf_concate(lh_surf_bold_files, LH_FS6_N_VERTEX)

			# Load volume BOLD data
			vol_bold_files = sorted(glob(os.path.join(BASE_DATA_DIR, sub, '*', 'Preprocess', 'vol', '*MNI2mm*.nii.gz')))
			vol_bold = get_vol_concate(vol_bold_files)

			# Compute FC maps for each seed
			for seed_name in SEED_LIST:
				seed_mask_file = os.path.join(SEED_MASK_DIR, f'{seed_name}_MNI2mm.nii.gz')
				seed_mask = ants.image_read(seed_mask_file)
				seed_mask_np = seed_mask.numpy()

				seed_vector = vol_bold[seed_mask_np == 1, :].mean(axis=0)

				# Right hemisphere
				fc_map = compute_fcmap(seed_vector, rh_surf_bold)
				fc_map = np.arctanh(fc_map)  # Fisher-Z transformation
				save_mgh(fc_map, f'{flag_save_path}/rh_Net_fcmap_{sub}_{seed_name}.mgh')

				# Left hemisphere
				fc_map = compute_fcmap(seed_vector, lh_surf_bold)
				fc_map = np.arctanh(fc_map)  # Fisher-Z transformation
				save_mgh(fc_map, f'{flag_save_path}/lh_Net_fcmap_{sub}_{seed_name}.mgh')
