# -*- coding: utf-8 -*-
# @Author  : zhangwei
# @Time    : 2024/3/26 AM 9:42
# @Func    : Fig1_bar

from utils.utils import *

if __name__ == '__main__':
    set_environ()

    """
    1. e.g. HC60 dataset
    """
    save_dir = '/path/to/HC60_single_roi_0809_stats'
    os.makedirs(save_dir, exist_ok=True)

    seed_list = ['Pu', 'GPe', 'GPi', 'STN', 'SN', 'VIM']
    plot_col = ['SCAN', 'Effector']
    stats_data = []
    for data_flag in ['HC']:
        for seed in seed_list:
            info = pd.read_csv(f'/path/to/Network_FC/BigPD_PD166_HC60_indiparc/{data_flag}_{seed}.csv')
            plot_info = info[plot_col]
            # save
            plot_info.to_csv(f"{save_dir}/{seed}_{data_flag}.csv", index=False)
            # stats
            SCAN_mean, SCAN_std = np.mean(plot_info['SCAN']), np.std(plot_info['SCAN'])
            Effector_mean, Effector_std = np.mean(plot_info['Effector']), np.std(plot_info['Effector'])
            r, p = calculate_ttest_paired(plot_info['SCAN'], plot_info['Effector'])
            stats_data.append([seed, f"{SCAN_mean:.3f}, {SCAN_std:.3f}", f"{Effector_mean:.3f}, {Effector_std:.3f}", f"{r:.3f}", f"{p:.4f}"])
            # plot
            plot_info = pd.melt(plot_info, var_name='Net', value_name="value")
            bar_plot(plot_info, 'Net', 'value', -0.1, 0.2, (1.5, 2), 0.1, save_dir, f"{seed}_{data_flag}")
        # save stats info
        stats_info = pd.DataFrame(data=stats_data, columns=['ROI', 'SCAN_FC(mean,std)', 'Effector_FC(mean,std)', 'r_value', 'P_value'])
        stats_info.to_csv(f"{save_dir}/{data_flag}_stats.csv", index=False)

    """
    3. e.g. PD65 dataset
    """
    save_dir = '/path/to/BigPD65_single_roi_0809_stats'
    os.makedirs(save_dir, exist_ok=True)

    seed_list = ['Pu', 'GPe', 'GPi', 'STN', 'SN', 'VIM']
    plot_col = ['SCAN', 'Effector']
    stats_data = []
    for data_flag in ['PD']:
        for seed in seed_list:
            info = pd.read_csv(f'/path/to/Network_FC/BigPD_PD65_HC60_indiparc/{data_flag}_{seed}.csv')
            info['Effector'] = (info['Foot-SM'] + info['Hand-SM'] + info['Face-SM']) / 3
            plot_info = info[plot_col]
            # save
            plot_info.to_csv(f"{save_dir}/{seed}_{data_flag}.csv", index=False)
            # stats
            SCAN_mean, SCAN_std = np.mean(plot_info['SCAN']), np.std(plot_info['SCAN'])
            Effector_mean, Effector_std = np.mean(plot_info['Effector']), np.std(plot_info['Effector'])
            r, p = calculate_ttest_paired(plot_info['SCAN'], plot_info['Effector'])
            stats_data.append([seed, f"{SCAN_mean:.3f}, {SCAN_std:.3f}", f"{Effector_mean:.3f}, {Effector_std:.3f}", f"{r:.3f}", f"{p:.4f}"])
            # plot
            plot_info = pd.melt(plot_info, var_name='Net', value_name="value")
            bar_plot(plot_info, 'Net', 'value', -0.1, 0.2, (1.5, 2), 0.1, save_dir, f"{seed}_{data_flag}")
        # save stats info
        stats_info = pd.DataFrame(data=stats_data, columns=['ROI', 'SCAN_FC(mean,std)', 'Effector_FC(mean,std)', 'r_value', 'P_value'])
        stats_info.to_csv(f"{save_dir}/{data_flag}_stats.csv", index=False)