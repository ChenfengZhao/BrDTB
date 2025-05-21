from pyteomics import mass
import numpy as np
from matchms.filtering import normalize_intensities
from matchms import Spectrum
from matchms.similarity import CosineGreedy
import pandas as pd
import sys
import math
import os
import glob
from matchms import set_matchms_logger_level
from collections import defaultdict
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
from matchms.importing import load_from_mzxml
import re
import numba
from typing import Tuple

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# TODO MS1 +-0.5min; If the first two theorical peak match, add up the two matched experimental peak intensity. Rank it and pick the 3 scans/spectra with highest sumed intensity and then calculate similarity.

@numba.njit
def find_matches_custom(spec1_mz: np.ndarray, spec2_mz: np.ndarray, spec2_inten: np.ndarray, tolerance: float, shift: float = 0) -> Tuple[float, float]:
    """Faster search for matching peaks.
    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from
    low to high m/z). Return the total intensity of matched experimental peaks and matched peak number.

    Derived from matchms.similarity.spectrum_similarity_functions.find_matches. If there are more than 1 matched spec2 peaks are detected for the same spec1, pick the spec2 with highest intensity.

    Parameters
    ----------
    spec1_mz:
        Theoretical Spectrum peak m/z values as numpy array. Peak mz values must be ordered.
    spec2_mz:
        Experimental Spectrum peak m/z values as numpy array. Peak mz values must be ordered.
    spec2_inten:
        Experimental 
    tolerance
        Peaks will be considered a match when <= tolerance appart.
    shift
        Shift peaks of second spectra by shift. The default is 0.

    Returns
    -------
    tot_inten:
        Total intensity of experimental peaks
    matched_num:
        The number of matched peaks

    """
    lowest_idx = 0
    tot_inten = 0
    matched_num = 0
    # matches = []
    for peak1_idx in range(spec1_mz.shape[0]):
        mz = spec1_mz[peak1_idx]
        low_bound = mz - tolerance
        high_bound = mz + tolerance
        spec2_inten_list = []
        for peak2_idx in range(lowest_idx, spec2_mz.shape[0]):
            mz2 = spec2_mz[peak2_idx] + shift
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                lowest_idx = peak2_idx + 1
            else:
                # matches.append((peak1_idx, peak2_idx))
                # matches.append((spec1_mz[peak1_idx], spec2_mz[peak2_idx]))
                # tot_inten += spec2_inten[peak2_idx]
                spec2_inten_list.append(spec2_inten[peak2_idx])
                # matched_num += 1
        
        if spec2_inten_list:
            # tot_inten += max(spec2_inten_list)
            tot_inten += sum(spec2_inten_list)
            matched_num += 1
    
    # print("matches:", matches)
    return tot_inten, matched_num

def creat_path(path):
    """Creat path+folder if not exist. Do nothing if path exists

    Parameters
    ----------
    path : str
        path + folder_name
    """
    isExists = os.path.exists(path)

    if not isExists:
        os.makedirs(path)
        # print("path generated:", path)
    else:
        # print("path exists:", path)
        pass

def find_filter_peaks(args, mz_lst, intensity_lst, filter=False, sigma=1.0, min_height=0.0, threshold=0.0, debug=False):
    """Find the peaks in mz-intensity

    Parameters
    ----------
    args : Class
        Input arguments
    mz_lst : list
        mz list
    intensity_lst : list
        intensity list
    filter : bool, optional
        apply gaussian filter at the first begining if set True, by default False
    sigma : float, optional
        sigma used in Gaussian filter, by default 1.0
    min_height : float, optional
        minimum peaks' intensity, by default 0
    threshold : float, optional
        portion of minimum peaks' intensity to the maximum peaks' intensity, by default 0. Range: 0-1
    debug : bool, optional
        print important messages if set true, by default False.

    Return
    ------
    peak_mz_lst : list
        mz list of peaks
    peak_intensity_lst : list
        intensity list of peaks

    """

    # Apply a Gaussian filter with sigma if filter is set True.
    if filter:
        intensity_arr = gaussian_filter(intensity_lst, sigma=sigma)
    else:
        intensity_arr = np.array(intensity_lst)

    # find the max intensity and min intensity used to filter some peaks
    max_intensity = intensity_arr.max()
    min_intensity = max(min_height, threshold * max_intensity)

    if debug:
        print("INFO: max_intensity", max_intensity, "\n", "min_intensity:", min_intensity)
    
    # find all the peaks
    peak_idx_arr, find_peak_prop_dict = find_peaks(intensity_arr)

    if debug:
        print("INFO: indices of all the peaks::", peak_idx_arr)
        print("INFO: find_peak_prop_dict:", find_peak_prop_dict)
    

    peak_mz_lst = [mz_lst[i] for i in peak_idx_arr]
    peak_intensity_lst = [intensity_arr[i] for i in peak_idx_arr]

    return peak_mz_lst, peak_intensity_lst

def save_vlineFig2Pdf(pdf, x, y, title, color = 'blue', xmin = None, xmax = None):
    """Used alone with "with PdfPages("two plots.pdf") as pdf" to save a figure to pdf

    Parameters
    ----------
    pdf : _type_
        pdf handler created by PdfPages()
    x : list or numpy.array
        x axis values
    y : list or numpy.array
        y axis values
    color : str
        color of the points and vlines
    title : str
        title of the figure
    """

    # first figure
    plt.figure()
    plt.scatter(x, y, color=color)  # Plot the points

    # Draw vertical lines from x-axis to each (x, y) point
    for xi, yi in zip(x, y):
        plt.vlines(xi, 0, yi, colors=color)

    if xmin:
        # print("setting xmin")
        plt.xlim(left=xmin)
    if xmax:
        # print("setting xmax")
        plt.xlim(right=xmax)
    
    plt.ylim(bottom=0)

    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.title(title)
    plt.grid(True)
    # plt.savefig('plot_with_lines.pdf') 
    pdf.savefig()
    plt.close()

def read_rst_summary_excel(fp):
    """read the result summary excel file (fp: file name + path). If fp doesn't exist, create a new excel file. Return the df of the result summary excel file

    Parameters
    ----------
    fp : str
        file name + path of the result summary excel file
    
    Return
    -------
    df : pd.DataFrame
    """
    if os.path.exists(fp):
        df = pd.read_excel(fp)
    else:
        data = {'Peptide_Name': [], 
                'z': [], 
                'Data_File': [],
                'Scan_Number':[],
                'Similarity_Score(Br)': [],
                'Similarity_Score(nonBr)': []}
        df = pd.DataFrame(data)
        df.to_excel(fp, index=False, engine='openpyxl')
    
    return df


def filter_peaks(mz_inten_lst, tolerance=0.0088):
    """filter adjacent peaks (mz, intensity) whose mz is in range of tolerance

    Parameters
    ----------
    mz_inten_lst : list of tuple
        list of tuple (mz, intensity)
    tolerance : float, optional
        mz tolerance of adjacent peaks, by default 0.0088

    Returns
    -------
    result : list of tuple
        list of tuple (mz, intensity)
    """
    # Sort by m/z value
    sorted_lst = sorted(mz_inten_lst, key=lambda x: x[0])
    
    result = []
    current_group = []
    
    for peak in sorted_lst:
        if not current_group:
            current_group.append(peak)
        else:
            # Compare with the first peak in the current group
            if abs(peak[0] - current_group[0][0]) <= tolerance:
                current_group.append(peak)
            else:
                # Keep the peak with the highest intensity
                max_peak = max(current_group, key=lambda x: x[1])
                result.append(max_peak)
                current_group = [peak]
    
    # Process last group
    if current_group:
        max_peak = max(current_group, key=lambda x: x[1])
        result.append(max_peak)
    
    return result

if __name__ == "__main__":

    set_matchms_logger_level("ERROR")

    # excel_file = "../data/20240802_BBP_Peptide_1_BrDTB.raw_20240819_Byonic.xlsx"
    excel_file = "./data/test_peptide_BrDTB.xlsx"
    data_path = "./data"
    delta_t = 0.5
    pep_mod_dict = {"+658.259" : "C30H43BrN8O4"} # peptide modification dict {modification mass: modification composition}
    abd_threshold = 0.005 # filter out the theoritical peaks that are lower than abd_threshold, reletive to the first theoritical peaks
    rst_path = "./results"
    rst_summary_fn = "result_summary.xlsx"

    # check whether the results summary excel file exists
    creat_path(rst_path)
    pep_done_df = read_rst_summary_excel(rst_path + "/" + rst_summary_fn)
    # get all the peptides that have been processed, and save it into a set
    pep_done_set = set()
    for index, row in pep_done_df.iterrows():
        pep_n = row['Peptide_Name']
        pep_z = row['z']
        pep_dn = row['Data_File']
        pep_done_set.add(pep_n + '-' + str(pep_z) + '-' + pep_dn)

    # information of target peptides 
    pep_info_df = pd.read_excel(excel_file, sheet_name="Spectra") 
    # print("pep_info_df:\n", pep_info_df)

    # Convert the excel_file to the following dicts 
    pep_rt_dict = defaultdict(list) # {peptide_charge: list of scan time}
    pep_fn_dict = defaultdict(str) # {peptide_charge: peptide file name}

    for index, row in pep_info_df.iterrows():
        # print(row)
        # print(row['Peptide\n< ProteinMetrics Confidential >'])
        pep_n = row['Peptide\n< ProteinMetrics Confidential >'] # peptide name
        pep_z = row['z'] # peptide charge
        pep_rt = row['Scan Time'] # peptide scan/retention time (ms2)
        pep_dn = row['Comment'].split('.')[0] # data file name

        pep_rt_dict[pep_n + '-' + str(pep_z) + '-' + pep_dn].append(pep_rt)
        # pep_fn_dict[pep_n + '-' + str(pep_z)] = row['Comment'].split('.')[0]
    
    # print("pep_rt_dict:", pep_rt_dict)
    # print("pep_fn_dict:",pep_fn_dict)

    # Find the data path + data name of all the mzXML files under data_path, including all the subpaths
    data_pn_lst = glob.glob(data_path+"/**/*.mzXML", recursive=True) # data path + data name of all the data files, including all the subpaths


    # Process each peptide with specific charge
    for pep_n_z_dn, tgt_t_lst in pep_rt_dict.items():

        pep_n = pep_n_z_dn.split('-')[0]
        pep_z = int(pep_n_z_dn.split('-')[1])
        tgt_fn = pep_n_z_dn.split('-')[2]

        # if pep_n_z_dn exists in pep_done_set, skip this pep
        if pep_n_z_dn in pep_done_set:
            print("INFO: Peptide: %s, charge: %s, data File: %s has already been processed. Skip this peptide!" % (pep_n, str(pep_z), tgt_fn))
            continue
        # Create a new row for result summary excel file
        row_done = {'Peptide_Name': pep_n, 'z': str(pep_z), 'Data_File': tgt_fn, 'Scan_Number':[], 'Similarity_Score(Br)': [], 'Similarity_Score(nonBr)': []}


        print("INFO: Processing peptide: %s, charge: %s, data File: %s......"% (pep_n, str(pep_z), tgt_fn))
        # Find the matched data file which has the same name as peptide name, regardless of big/small capital
        matched_file_lst = [fpn for fpn in data_pn_lst if os.path.basename(fpn).split(".")[0].lower() == tgt_fn.lower()]

        if len(matched_file_lst) != 1:
            # print("matched_file_lst:", matched_file_lst)
            if not matched_file_lst:
                print("ERROR: no matched data file found!")
            else:
                print("ERROR: More than 2 matched data files found!", matched_file_lst)
            exit()

        # print("matched_file_lst", matched_file_lst)
        matched_fpn = matched_file_lst[0]
        print("INFO: found matched file!", matched_fpn)

        # Read ms1 scans of the matched data file
        ms1_spectra_lst = list(load_from_mzxml(matched_fpn, ms_level=1))

        # Calcualte target scan/retention time tgt_t
        # tgt_t_lst = pep_rt_dict[pep_n_z]
        if not tgt_t_lst:
            print("ERROR: scan time list of %s in charge %s is empty!" % (pep_n, str(pep_z)))
            exit()
        # Calcuate the averge of tgt_t_lst as tgt_t of ms2
        tgt_t = sum(tgt_t_lst)/len(tgt_t_lst)

        # Get the peptide composition from pep_n
        # pep_n = pep_n_z.split('-')[0] # raw peptide name in excel
        # pep_z = int(pep_n_z.split('-')[1]) # peptide charge in excel
        pep_match = re.search(r'\.([A-Z]+)\[(.*?)\]([A-Z]+)\.', pep_n)
        if not pep_match:
            print("ERROR: " + pep_n + " doesn't match the pattern!")
            exit()
        pep_seq = pep_match.group(1) + pep_match.group(3)
        pep_mod = pep_match.group(2)
        pep_mod_comp = pep_mod_dict[pep_mod]
        # print('pep_seq:', pep_seq)
        # print('pep_mod_comp:', pep_mod_comp)
        pep_comp = mass.Composition(sequence=pep_seq) + mass.Composition(formula=pep_mod_comp)
        # print('pep_comp:', pep_comp)

        # Get theoritical isotopologue peaks: list of tuple [(mz, abundance)]
        # @ToDo: nomalization to the first peak and then overall threshold = 0.05
        # @ToDo: mass -> m/z: (mass + 3 * H+)/3
        mz_abd_lst = [(mass.calculate_mass(composition=isotopologue, charge=pep_z), abundance) for isotopologue, abundance in mass.isotopologues(pep_comp, report_abundance=True, overall_threshold=1e-5)]
        mz_abd_lst.sort()
        # print('mz_abd_lst:', mz_abd_lst)

        # nomalize abundance based on the first abundance
        if len(mz_abd_lst):
            abd_base = mz_abd_lst[0][1]
            mz_abd_norm = [(mz, abundance/abd_base) for mz, abundance in mz_abd_lst if abundance/abd_base >= abd_threshold]
            # print("mz_abd_norm:", mz_abd_norm)
        else:
            print("ERROR: No isotopologue found!", file=sys.stderr)
            sys.exit(1)

        # filter peaks whose mz is too close to neighbors
        mz_abd_norm = filter_peaks(mz_abd_norm)

        mz_lst, abd_norm_lst = np.array(mz_abd_norm).T
        # print('mz_lst:', mz_lst, type(mz_lst))
        # print('abd_norm_lst:', abd_norm_lst, type(abd_norm_lst))

        mz_max = mz_lst.max()
        mz_min = mz_lst.min()
        # calculated Spectrum objects
        calc_spectrum = Spectrum(mz=mz_lst, intensities=abd_norm_lst, metadata={"id": "calculated"})
        # normalize the theoretical spectra before comparison
        calc_spectrum = normalize_intensities(calc_spectrum)


        # Search ms1 spectra to find the sepctra whose retention time is within [tgt_t - 0.5, tgt_t + 0.5]. Calculate the similarity of all the selected spectra, and pick 3 spectra with highest similarity score.
        score_spec_lst = [] # [(score, exp_spectrum)]
        for spectrum in ms1_spectra_lst:
            # filter out spectra whose retention time is not within [tgt_t - 0.5, tgt_t + 0.5]
            spt_rt = spectrum.metadata['retention_time']
            if (spt_rt < tgt_t - delta_t) or (spt_rt > tgt_t + delta_t):
                continue

            # filter out the peaks that is too far way from theoritical peaks
            spt_mz = spectrum.mz
            spt_inten = spectrum.intensities
            mask = (spt_mz >= (mz_min - 0.07)) & (spt_mz <= (mz_max + 0.07))
            filt_spt_mz = spt_mz[mask]
            filt_spt_inten = spt_inten[mask]

            # create experiment spectrum
            exp_spectrum = Spectrum(mz=filt_spt_mz, intensities=filt_spt_inten, metadata=spectrum.metadata)

            # normalize the experimental spectra before comparison
            exp_spectrum = normalize_intensities(exp_spectrum)

            # Compute Similarity
            similarity_measure = CosineGreedy(tolerance=0.07)  # Tolerance for matching peaks
            score = similarity_measure.pair(calc_spectrum, exp_spectrum)

            # score_spec_lst
            # print("score['score']", score['score'], type(score['score']))
            # score_spec_lst.append((score['score'], spectrum.metadata['scan_number']))
            score_spec_lst.append((score['score'], exp_spectrum))

        score_spec_lst.sort(reverse=True, key=lambda x: x[0])
        # print("score_spec_lst:\n", score_spec_lst)
        score_spec_lst = score_spec_lst[:3]


        # ###### Calculate similarity score of non-Br Composition ######
        # calculate the theorical spectrum of non-Br Composition
        nonBr_pep_comp = pep_comp - mass.Composition(formula='Br') + mass.Composition(formula="C6H7")
        print("INFO: nonBr peptide composition:", nonBr_pep_comp)

        nonBr_mz_abd_lst = [(mass.calculate_mass(composition=isotopologue, charge=pep_z), abundance) for isotopologue, abundance in mass.isotopologues(nonBr_pep_comp, report_abundance=True, overall_threshold=1e-5)]
        nonBr_mz_abd_lst.sort()

        # nomalize nonBr abundance based on the first abundance
        if len(nonBr_mz_abd_lst):
            abd_base = nonBr_mz_abd_lst[0][1]
            nonBr_mz_abd_norm = [(mz, abundance/abd_base) for mz, abundance in nonBr_mz_abd_lst if abundance/abd_base >= abd_threshold]
            # print("nonBr_mz_abd_norm:", nonBr_mz_abd_norm)
        else:
            print("ERROR: No isotopologue found!", file=sys.stderr)
            sys.exit(1)

        # filter peaks whose mz is too close to neighbors
        nonBr_mz_abd_norm = filter_peaks(nonBr_mz_abd_norm)
        
        nonBr_mz_lst, nonBr_abd_norm_lst = np.array(nonBr_mz_abd_norm).T

        # calculated Spectrum objects
        nonBr_calc_spectrum = Spectrum(mz=nonBr_mz_lst, intensities=nonBr_abd_norm_lst, metadata={"id": "calculated"})
        # normalize the theoretical spectra before comparison
        nonBr_calc_spectrum = normalize_intensities(nonBr_calc_spectrum)

        # calculate the similarity score between the theorical spectrum of non-Br compisition and the 3 picked experimental spectra
        nonBr_score_lst = [] # scores of nonBr theoretical
        for _, exp_spectrum in score_spec_lst:
            # Compute Similarity
            similarity_measure = CosineGreedy(tolerance=0.07)  # Tolerance for matching peaks
            nonBr_score = similarity_measure.pair(nonBr_calc_spectrum, exp_spectrum)
            # print("nonBr_score:", nonBr_score)
            nonBr_score_lst.append(nonBr_score)
            
        
        # save figures to 3 pdf files. each pdf file contains an experimental spectrum, Br theoretical spectrum, and a non-Br theoretical specturm.
        pdf_base_n = pep_seq + '-' + str(pep_z) + '-' + pep_dn
        xmin = calc_spectrum.mz[0] - 0.1
        xmax = calc_spectrum.mz[-1] + 0.1
        for i, (Br_score, exp_spectrum) in enumerate(score_spec_lst):
            scan_num = exp_spectrum.metadata["scan_number"]
            pdf_n = pdf_base_n + '-' + scan_num + '.pdf'
            nonBr_score = nonBr_score_lst[i]['score']
            with PdfPages(rst_path + "/" + pdf_n) as pdf:

                # experimental figure
                save_vlineFig2Pdf(pdf, exp_spectrum.mz, exp_spectrum.intensities, "Normalized Experimental Spectrum, Scan Number: " + scan_num, color="blue", xmin=xmin, xmax=xmax)

                # Br theoretical figure
                save_vlineFig2Pdf(pdf, calc_spectrum.mz, calc_spectrum.intensities, "Normalized Theoretical Spectrum (Br), Similarity: %.3f" % Br_score, color="blue", xmin=xmin, xmax=xmax)

                # Non-Br theoretical figure
                save_vlineFig2Pdf(pdf, nonBr_calc_spectrum.mz, nonBr_calc_spectrum.intensities, "Normalized Theoretical Spectrum (w/o Br), Similarity: %.3f" % nonBr_score, color="blue", xmin=xmin, xmax=xmax)

                print("INFO: scan number: " + scan_num + ", Br similarity score: %.3f, nonBr similarity score:  %.3f" %(Br_score, nonBr_score))
            
            row_done['Scan_Number'].append(scan_num)
            row_done['Similarity_Score(Br)'].append(Br_score)
            row_done['Similarity_Score(nonBr)'].append(nonBr_score)
        
        row_done['Scan_Number'] = ",".join(map(str, row_done['Scan_Number']))
        row_done['Similarity_Score(Br)'] = ",".join([f"{x:.3f}" for x in row_done['Similarity_Score(Br)']])
        row_done['Similarity_Score(nonBr)'] = ",".join([f"{x:.3f}" for x in row_done['Similarity_Score(nonBr)']])
        pep_done_df = pd.read_excel(rst_path + "/" + rst_summary_fn)
        pep_done_df = pd.concat([pep_done_df, pd.DataFrame([row_done])], ignore_index=True)
        pep_done_df.to_excel(rst_path + "/" + rst_summary_fn, index=False, engine='openpyxl')




                



