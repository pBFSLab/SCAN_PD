# -*- coding: utf-8 -*-
# @Author  : zhangwei
# @Time    : 2023/6/25 PM 4:01
# @Func    : S2_calcu_network_FC

from utils.utils import *

if __name__ == '__main__':
	"""
	1. Cortex mask
	"""
	parc_labels = ['10', '11', '17', '18']
	Net_names = ['Hand-SM', 'Face-SM', 'Foot-SM', 'SCAN']

	HC_lh_parc_base_path = '/path/to/HC_lh_parc_base_path/{}/*/Iter_10_sm4/lh.parc_result.annot'
	HC_rh_parc_base_path = '/path/to/HC_rh_parc_base_path/{}/*/Iter_10_sm4/rh.parc_result.annot'
	PD_lh_parc_base_path = '/path/to/PD_lh_parc_base_path/{}/*/Iter_10_sm4/lh.parc_result.annot'
	PD_rh_parc_base_path = '/path/to/PD_rh_parc_base_path/{}/*/Iter_10_sm4/rh.parc_result.annot'

	"""
	2. FC data
	"""
	PD_control_info = pd.read_csv(
		'/path/to/PD_control_info.csv')  # 60
	PD_control_data_list = PD_control_info['Subject_ID'].tolist()

	# PD patient
	PD_patient_info = pd.read_csv(
		'/path/to/PD_patient_info.csv')  # 65
	PD_patient_data_list = PD_patient_info['Subject_ID'].tolist()
	print(f"HC:{len(PD_control_data_list)}, PD:{len(PD_patient_data_list)}")

	"""
	3. Calculate Network FC
	"""
	workdir = f'/path/to/workdir/Network_FC/BigPD_PD{len(PD_patient_data_list)}_HC{len(PD_control_data_list)}_indiparc'
	if not os.path.exists(workdir):
		os.makedirs(workdir)

	# (1) Save all ROI and Cortex FC mat files
	seed_list = ['Pu', 'GPe', 'GPi', 'STN', 'SN', 'VIM']

	for seed in seed_list:
		print(f"seed:{seed}")
		data_list = []
		HC_lh_base_path = '/path/to/HC_lh_base_path/lh_Net_fcmap_{}_{}.mgh'
		HC_rh_base_path = '/path/to/HC_rh_base_path/rh_Net_fcmap_{}_{}.mgh'
		PD_lh_base_path = '/path/to/PD_lh_base_path/lh_Net_fcmap_{}_{}.mgh'
		PD_rh_base_path = '/path/to/PD_rh_base_path/rh_Net_fcmap_{}_{}.mgh'

		"""
		HC
		"""
		HC_grp_lh = np.zeros((40962, len(PD_control_data_list)))
		HC_grp_rh = np.zeros((40962, len(PD_control_data_list)))
		grp_num = 0
		for idx, sub in enumerate(PD_control_data_list):
			lh_temp = HC_lh_base_path.format(sub, seed)
			rh_temp = HC_rh_base_path.format(sub, seed)
			fc_map_lh = nib.load(lh_temp).get_fdata().flatten()
			fc_map_rh = nib.load(rh_temp).get_fdata().flatten()

			HC_grp_lh[:, idx] = fc_map_lh
			HC_grp_rh[:, idx] = fc_map_rh

		# Save all FC data
		save_array_to_mat_v2(HC_grp_lh, 'grp_lh', HC_grp_rh, 'grp_rh', f'{workdir}/{seed}_HC_grp_data.mat')

		"""
		PD
		"""
		PD_grp_lh = np.zeros((40962, len(PD_patient_data_list)))
		PD_grp_rh = np.zeros((40962, len(PD_patient_data_list)))
		grp_num = 0
		for idx, sub in enumerate(PD_patient_data_list):
			lh_temp = PD_lh_base_path.format(sub, seed)
			rh_temp = PD_rh_base_path.format(sub, seed)

			fc_map_lh = nib.load(lh_temp).get_fdata().flatten()
			fc_map_rh = nib.load(rh_temp).get_fdata().flatten()

			PD_grp_lh[:, idx] = fc_map_lh
			PD_grp_rh[:, idx] = fc_map_rh

		save_array_to_mat_v2(PD_grp_lh, 'grp_lh', PD_grp_rh, 'grp_rh', f'{workdir}/{seed}_PD_grp_data.mat')

	# (3) Calculate FC intensity in Gorden surface17

	for seed in seed_list:
		"""
		FC Calculate, indiparc
		"""
		print(seed)
		HC_grp_lh, HC_grp_rh = load_mat_to_array_v2(f'{workdir}/{seed}_HC_grp_data.mat', 'grp_lh', 'grp_rh')
		PD_grp_lh, PD_grp_rh = load_mat_to_array_v2(f'{workdir}/{seed}_PD_grp_data.mat', 'grp_lh', 'grp_rh')

		HC_fc_dict = {}
		PD_fc_dict = {}
		for i, parc in enumerate(parc_labels):
			net = Net_names[i]
			parc = int(parc)
			print(f'net:{net}')

			# PD_patient
			PD_fc = []
			for idx, sub in enumerate(PD_patient_data_list):
				lh_annot_p = glob(PD_lh_parc_base_path.format(sub))[0]
				rh_annot_p = glob(PD_rh_parc_base_path.format(sub))[0]

				lh_labels, _, _ = nib.freesurfer.read_annot(lh_annot_p)
				rh_labels, _, _ = nib.freesurfer.read_annot(rh_annot_p)
				lh_mgh = np.zeros_like(lh_labels)
				rh_mgh = np.zeros_like(rh_labels)
				lh_mgh[lh_labels == parc] = 1
				rh_mgh[rh_labels == parc] = 1
				# Calculate FC intensity
				sub_fc = cortex_roi_fc_mean(PD_grp_lh[:, idx], PD_grp_rh[:, idx], lh_mgh, rh_mgh)
				PD_fc.append(sub_fc)

			# HC
			HC_fc = []
			for idx, sub in enumerate(PD_control_data_list):
				lh_annot_p = glob(HC_lh_parc_base_path.format(sub))[0]
				rh_annot_p = glob(HC_rh_parc_base_path.format(sub))[0]

				lh_labels, _, _ = nib.freesurfer.read_annot(lh_annot_p)
				rh_labels, _, _ = nib.freesurfer.read_annot(rh_annot_p)
				lh_mgh = np.zeros_like(lh_labels)
				rh_mgh = np.zeros_like(rh_labels)
				lh_mgh[lh_labels == parc] = 1
				rh_mgh[rh_labels == parc] = 1

				# Calculate FC intensity
				sub_fc = cortex_roi_fc_mean(HC_grp_lh[:, idx], HC_grp_rh[:, idx], lh_mgh, rh_mgh)
				HC_fc.append(sub_fc)

			PD_fc_dict[net] = PD_fc
			HC_fc_dict[net] = HC_fc

		HC_fc_df = pd.DataFrame(HC_fc_dict)
		PD_fc_df = pd.DataFrame(PD_fc_dict)

		HC_fc_df.to_csv(f'{workdir}/HC_{seed}.csv', index=False)
		PD_fc_df.to_csv(f'{workdir}/PD_{seed}.csv', index=False)
