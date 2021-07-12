import pyteomics.mzxml
import pandas as pd
import numpy as np
import tqdm
import scipy.signal
from sklearn.linear_model import LinearRegression

def check_local_max(spectrum, df_clip, min_int, max_int, around=2):
    
    bo_ar = (df_clip.iloc[spectrum]['m/z array'] > min_int) & (df_clip.iloc[spectrum]['m/z array'] < max_int)
    bo_ar2 = (((df_clip.iloc[spectrum]['m/z array'] > min_int-around) & (df_clip.iloc[spectrum]['m/z array'] < min_int)) |
              ((df_clip.iloc[spectrum]['m/z array'] > max_int) & (df_clip.iloc[spectrum]['m/z array'] < max_int+around)))    
    
    if df_clip.iloc[spectrum]['intensity array'][bo_ar].max() > df_clip.iloc[spectrum]['intensity array'][bo_ar2].max():
        return True
    else:
        return False

def find_peaks(chromats, peak_ids=['C23', 'C25', 'C27', 'C29'], start_time=12,
               peak_data={'C23':{'min_ion': 324, 'max_ion': 325, 'min_time':13, 'max_time':17, 'prom':250, 'width': 3},
                          'C25':{'min_ion': 352, 'max_ion': 353, 'min_time':15, 'max_time':19, 'prom':250, 'width': 3},
                          'C27':{'min_ion': 380, 'max_ion': 381, 'min_time':17, 'max_time':21, 'prom':250, 'width': 3},
                          'C29':{'min_ion': 408, 'max_ion': 409, 'min_time':19, 'max_time':23, 'prom':250, 'width': 3}},
              min_prom=230, prom_stride=10):
    
    """
    Find location of compounds of interest based on GCMS data and its diagnostic ion.
    Parameters
    ----------
    chromats : array like
        List or array like object with the absolute path to a set of GCMS files (mzXML format) to find
        peaks of interst in
    peak_ids : array like
        List or array like object with the compound ID corresponding to a dictionary of
        info on the diagnostic ion to use to locate that compound
    start_time : float
        What time in the chromatogram to start looking for compounds. Can speed
        up processing to ignore lections of the GCMS data.
    peak_data : dictionary
        A dictionary of dictionaries whose keys are peak_ids (a nickname for a given compound of interest).
        Each entry of the dictionary is a dictionary with info on how to ID the compound. The required
        entries of each dictionary are 'min_ion', the minimum M/Z value to consider as the diagnostic ion
        for the compound; 'max_ion', the maximal M/Z value to consider as the diagnostic ion
        for the compound; 'min_time', the earliest retention time to start looking for the diagnostic ion;
        'max_time', the maximal retetion time to look for the ion; 'prom', the starting peak prominance to use
        when looking for the peak in the diagnostic ion; 'width', the width (number of points) required for a
        peak call to be made.
    min_prom : float
        If the peak is not found for a given compound at a given peak prominance, the algorithm decreases
        the required prominance to make a call down to the min_prom defined here.
    prom_stride : 
        Amount to decrease the required peak prominance for a given compound call by before attempting again
        to find a peak.
    Returns
    -------
    peaks : pd.DataFrame
        DataFrame with retention time, peak intensity, peak_id, GCMS file, left peak base location, and right peak
        base location. If a peak is not found for a given compound, the values are zero for all features of the
        peak.
    """
    
    
    peaks = []
    for chromat in tqdm.tqdm(chromats):
    
        df_chrom = pd.DataFrame(pyteomics.mzxml.read(chromat))
        
        df_chrom = df_chrom.loc[df_chrom['retentionTime'] > start_time]

        for peak_id in peak_ids:
            min_int, max_int = (peak_data[peak_id]['min_ion'], peak_data[peak_id]['max_ion'])
            prom, wid = (peak_data[peak_id]['prom'], peak_data[peak_id]['width'])
            t1, t2 = (peak_data[peak_id]['min_time'], peak_data[peak_id]['max_time'])
            df_clip = df_chrom.loc[(df_chrom['retentionTime']>t1) & (df_chrom['retentionTime']<t2)]

            intensity = []
            for i in range(len(df_clip)):
                bo_ar = (df_clip.iloc[i]['m/z array'] > min_int) & (df_clip.iloc[i]['m/z array'] < max_int)
                if bo_ar.sum() > 0:
                    intensity.append((df_clip.iloc[i]['retentionTime'], df_clip.iloc[i]['intensity array'][bo_ar][0], i))
                    
            need_break=False
            while(prom>min_prom):
                try:
                    peak_loc = scipy.signal.find_peaks(np.array(intensity)[:, 1], prominence=prom, width=wid)
                except:
                    #print(intensity)
                    need_break=True
                    break
                if(len(peak_loc[0]) == 0):
                    #print(peak_id, prom, end='\t')
                    prom-=prom_stride
                else:
                    break
            
            if(need_break):
                break
            
            if len(peak_loc[0]) > 0:
                append_peak = False
                #print(peak_loc)
                peak = peak_loc[0][0]
                left = peak_loc[1]['left_bases'][0]
                right = peak_loc[1]['right_bases'][0]
                for peak_ind, peakl in enumerate(peak_loc[0]):
                    if check_local_max(intensity[peakl][2], df_clip, min_int, max_int):
                        peak = peakl
                        left = peak - int(peak_loc[1]['widths'][peak_ind])
                        right = peak + int(peak_loc[1]['widths'][peak_ind])
                        append_peak = True
                        break
                if append_peak:
                    peaks.append(np.array((np.array(intensity)[:, 0][peak], 
                                           np.array(intensity)[:, 1][peak],
                                          peak_id, 
                                          chromat,
                                          np.array(intensity)[:, 0][left],
                                          np.array(intensity)[:, 0][right])))
                else:
                    peaks.append(np.array((0, 0, peak_id, chromat, 0, 0)))
            else:
                peaks.append(np.array((0, 0, peak_id, chromat, 0, 0)))
    return(pd.DataFrame(peaks, columns=('ret time', 'intensity', 'identity', 'data', 'left_base', 'right_base'), dtype=np.double))

def align_trace(df_CHCs, reference, compare):
    
    times = []
    for chc in df_CHCs['identity'].unique():
        if  np.any(np.isclose(df_CHCs.loc[(df_CHCs['identity'] == chc) & ((df_CHCs['data'] == compare) | (df_CHCs['data'] == reference)), 'ret time'], 0)):
            continue
        times.append(df_CHCs.loc[(df_CHCs['identity'] == chc) & ((df_CHCs['data'] == reference) | (df_CHCs['data'] == compare))]['ret time'].values)
    if len(times) == 1:
        return(times[0][0]/times[0][1], 0, 1, times)
    if len(times) == 0:
        print('No peaks for '+ str(df_CHCs['identity'].unique()) + ' in the file '+ compare)
        return('skip', False, False, False)
    
    x = np.array(times)[:, 1].reshape((-1, 1))
    y = np.array(times)[:, 0]

    model = LinearRegression()

    model.fit(x, y)
    r_sq = model.score(x, y)

    return(model.coef_[0], model.intercept_, r_sq, times)

def calculate_baseline(baseline_file = "C:/GCMSsolution/Data/Julian_Wagner/20191217_scepto_lio_luct_from_20190401/20191217_hexane_blank_1_12172019_02.mzXML", start_time=12, start_delay=0, num_extrapolate=50, window=300):
    df_chrom = pd.DataFrame(pyteomics.mzxml.read(baseline_file))
    df_chrom = df_chrom.loc[df_chrom['retentionTime']>start_time]
    
    bg_smooth = np.zeros(len(df_chrom['intensity array']))
    bg_smooth[start_delay::] = df_chrom['intensity array'][start_delay::].apply(np.sum).rolling(window=window, win_type='general_gaussian', center=True).mean(width=100, power=1)

    bg_smooth[0:window+start_delay] = bg_smooth[start_delay+window:start_delay+window+num_extrapolate].mean()
    bg_smooth[len(bg_smooth)-window:len(bg_smooth)] = bg_smooth[len(bg_smooth)-window-num_extrapolate:len(bg_smooth)-window].mean()
    
    return(bg_smooth)
    #return(df_chrom['intensity array'].apply(np.sum))
    

def align_traces(df_CHCs, reference, baseline = False, norm_time=13, start_time=12):
    
    """
    Align a set of GCMS chromatograms based on the output of the compund detection
    done by the 'find_peaks' function.
    Parameters
    ----------
    df_CHCs : pandas DataFrame
        Output from `find_peaks` with retention time, peak intensity, peak_id, GCMS file,
        left peak base location, and right peak base location.
    reference : string
        Absolute path to GCMS file (mzXML format) to use as the reference to align the other chromatograms
        to. Must be one of the files that the peaks were detected in in the df_CHCs 
        DataFrame.
    baseline : string
        Absolute path to a GCMS file (mzXML format) to use as a baseline for background subtraction
        to impove visualization of GCMS chromatograms. Is only used for vizualization purposes.
    norm_time : float
        Retention time after which to normalize the GCMS trace for visualization purposes. Use
        0 if you wish to normalize based on the absolute max of the trace.
    start_time : float
        Gets rid of any portion of the GCMS trace before the start time.
    Returns
    -------
    df_traces : pandas DataFrame
        DataFrame with all the GCMS chromatograms in a tidy format (absolute as well as relative
        intensity), as well as a shifted and stretched retetion time that gives an aligned
        set of chromatograms.
    df_regressions : pandas DataFrame
        A DataFrame with 'slope' and 'intercept' needed to align a given chromatogram
        to the designated reference. Multiplying the retention time by the 'slope' and
        adding the 'intercept' gives the alinged chromatogram.
    """
    
    for i, comparison in tqdm.tqdm(enumerate(df_CHCs['data'].unique())):
        
        if not comparison == reference:
            coef, inter, r_sq, times = align_trace(df_CHCs, reference, comparison)
        else:
            coef, inter, r_sq = (1, 0, 1)
        
        if i == 0:
            df_regressions = pd.DataFrame({'data': comparison,
                                           'slope': coef,
                                           'intercept': inter,
                                           'r_sq': r_sq}, index=[0])
        else:
            df_regressions = df_regressions.append(pd.DataFrame({'data': comparison,
                                           'slope': coef,
                                           'intercept': inter,
                                           'r_sq': r_sq}, index=[i]))
        
        #coef, inter, r_sq = (1, 0, 1)
        
        if coef == 'skip':
            continue
            
            
        df_chrom = pd.DataFrame(pyteomics.mzxml.read(comparison))
        df_chrom = df_chrom.loc[df_chrom['retentionTime']>start_time]
        
        
        if baseline and not baseline[comparison] == 'NA':
            base = calculate_baseline(baseline_file=baseline[comparison], start_time=start_time)
        else:
            base = np.zeros(len(df_chrom['intensity array']))
        
        df_chrom['intensity minus base'] = df_chrom['intensity array'].apply(np.sum) - base#df_base['intensity array'].apply(np.sum).rolling(20).median()
        df_chrom['intensity minus base'] = df_chrom['intensity minus base'] - df_chrom['intensity minus base'].min()
        df_chrom['intensity norm'] = df_chrom['intensity minus base']/df_chrom.loc[df_chrom['retentionTime']>norm_time, 'intensity minus base'].max()
        
        if i == 0:
            df_traces = pd.DataFrame({'ret time': df_chrom['retentionTime'],
                                      'chromat':df_chrom['intensity array'].apply(np.sum),
                                      'data': comparison,
                                      'slope': coef,
                                      'intercept': inter,
                                      'r_sq': r_sq,
                                      'chromat_norm': df_chrom['intensity norm']})
        else:
            df_traces = df_traces.append(pd.DataFrame({'ret time': df_chrom['retentionTime'],
                                                       'chromat':df_chrom['intensity array'].apply(np.sum),
                                                       'data': comparison,
                                                       'slope': coef,
                                                       'intercept': inter,
                                                       'r_sq': r_sq,
                                                       'chromat_norm': df_chrom['intensity norm']}))
            
    df_traces['allign ret time'] = df_traces['ret time']*df_traces['slope'] + df_traces['intercept']
    return(df_traces, df_regressions)