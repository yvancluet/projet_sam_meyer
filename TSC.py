import os
import sys
import configparser
import pandas as pd
import numpy as np
import collections as col
from shutil import copy

###########################################################
#                       Functions                         #
###########################################################

def create_config_file(config, config_file_path, TSS_file, SIGMA_0, RNAPS_NB):  # DELTA_X, D, J_0, SIGMA_0, RNAPS_NB,
    # Create the config file
    config.set('INPUTS', 'TSS', str(TSS_file))
    config.set('SIMULATION', 'SIGMA_0', str(SIGMA_0))
    config.set('SIMULATION', 'RNAPS_NB', str(RNAPS_NB))
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)


# Read the config files
def read_config_file(path):
    config = configparser.ConfigParser()
    # to preserve capital letters
    config.optionxform = str
    # Update INI file without removing comments
    config = configparser.ConfigParser(allow_no_value=True)
    if not os.path.exists(path):
        print("Input file was not found at path %s" % path)
        sys.exit(1)
    config.read(path)
    return config


###################### Reading files ######################

# you can combine those two functions
def load_gff(filename):
    gff_df_raw = pd.read_csv(filename, sep='\t', comment='#', header=0)
    return gff_df_raw


def load_tab_file(filename):
    data = pd.read_csv(filename, sep='\t', header=0)
    return data


def str2num(s):
    s[s == '+'] = 1  # True
    s[s == '-'] = -1  # False
    return s


def get_tr_nbr_csv(csv_file):
    csv_tr_nbr = pd.read_csv(csv_file, sep='\t', header=None)
    tr_nbr = csv_tr_nbr.values
    return tr_nbr.flatten()


######################## Others ###########################

# Get the genome size from the header of gff file (befor renaming it)
def get_genome_size(gff_df):
    genome_size = int(gff_df.columns[4]) - int(gff_df.columns[3])
    return genome_size


# Rename the header (columns names)
def rename_gff_cols(gff_df):
    names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df.columns = names
    return gff_df


# calculate the initiation rate
def f_init_rate(tr_prob, sig, sigma_t, epsilon, m):
    tr_prob_sig = tr_prob * np.exp((1 / (1 + np.exp((sig - sigma_t) / epsilon))) * m)
    return tr_prob_sig


def is_gene_crossing_origin(ts_sorted, genome_size):
    """
    Check if there's a gene crossing the origin, there are four
    cases depending on the orientations of the first and last gene:
        1. ___G++(Ori)++____G++++__
        2. ___G++(Ori)++____G----__
        3. ___G--(Ori)--____G++++__
        4. ___G--(Ori)--____G----__
    Args:
        :param ts_sorted: pandas DataFrame of sorted Transcripts (Start and End) Sites.
        :param genome_size: genome size

    Returns:
        gene_crossing_origin (boolean) : true or false
    """
    gene_crossing_origin = False
    last_gene_size = 0

    # getting genes strands
    first_line_strand = ts_sorted['TUorient_x'].values[0]
    second_line_strand = ts_sorted['TUorient_x'].values[1]
    # getting genes tss positions
    first_line_tss_pos = ts_sorted['TSS_pos'].values[0]
    second_line_tss_pos = ts_sorted['TSS_pos'].values[1]

    # getting genes tts positions
    first_line_tts_pos = ts_sorted['TTS_pos'].values[0]
    #second_line_tts_pos = ts_sorted['TTS_pos'].values[1]

    # case : ___G++(Ori)++____G++++__  OR  case : ___G++(Ori)++____G----__
    if (first_line_strand == '+' and second_line_strand == '+') or (first_line_strand == '+' and second_line_strand == '-'):
        # equivalent to: second_line_tss_pos > first_line_tts_pos and first_line_tss_pos > second_line_tss_pos
        gene_crossing_origin = first_line_tts_pos < second_line_tss_pos < first_line_tss_pos
        last_gene_size = genome_size - first_line_tss_pos + first_line_tts_pos

    # case : ___G--(Ori)--____G++++__ OR case : ___G--(Ori)--____G----__
    elif (first_line_strand == '-' and second_line_strand == '+') or (first_line_strand == '-' and second_line_strand == '-'):
        # equivalent to: second_line_tss_pos > first_line_tss_pos and second_line_tss_pos < first_line_tts_pos
        gene_crossing_origin = first_line_tss_pos < second_line_tss_pos < first_line_tts_pos
        last_gene_size = genome_size - first_line_tts_pos + first_line_tss_pos

    return gene_crossing_origin, last_gene_size


def sort_by(all_tss_pos, all_tts_pos):
    """
    Sort either by TSS positions or TTS ones
    depending on where is the minimum value
    Args:
        all_tss_pos (list): list of all TSS position extracted from the 'tss.dat' file.
        all_tts_pos (list): list of all TTS position extracted from the 'tts.dat' file.
    Returns:
        col_name (string) : which will be either "TSS_pos" or "TSS_pos"
    """
    min_tss = min(all_tss_pos)
    min_tts = min(all_tts_pos)
    col_name = ""
    if min_tss < min_tts:
        col_name = "TSS_pos"
    else:
        col_name = "TTS_pos"
    return col_name


def clean_sorted_ts(ts):
    """
    Clean the sorted TS dataframe.
    After merging tts and tts file generate all the possible
    combinations between the two, this function cleans the
    unnecessary rows and return the TS dataframe ready to be used
    """
    # keep only the rows
    ts = ts[
        # where TSS and TTS have the same orientation
        (ts["TUorient_x"] == ts["TUorient_y"])
        &  # and we get rid of the TSS/TSS and TTS/TTS regions
        (
                (ts["TSS_pos"] < ts["TTS_pos"]) & (ts["TUorient_x"] == "+")
                |  # or
                (ts["TSS_pos"] > ts["TTS_pos"]) & (ts["TUorient_x"] == "-")
        )
        ]
    # reset the index (after removing the wrong rows)
    ts.reset_index(inplace=True)

    # this loop is used to get rid of rows that
    # doesn't take into account the Poff value, example:
    # _S____T=1____TT=1___ => S ---> T [CORRECT!]
    #                      => S ---> TT [WRONG!]
    # 'S': Start Site, 'T=1': Termination Site with Poff = 1
    # 'S ---> T' and 'S ---> TT': possible transcripts
    for ix, row in ts.iterrows():
        current_tss_pos = ts['TSS_pos'][ix]
        privious_tss_pos = ts['TSS_pos'].shift(1)[ix]
        current_tu_orient = ts['TUorient_x'][ix]
        privious_tu_orient = ts['TUorient_x'].shift(1)[ix]
        privious_poff = ts['TTS_proba_off'].shift(1)[ix]
        if (
                ix > 0 and
                current_tss_pos == privious_tss_pos and
                current_tu_orient == privious_tu_orient and
                privious_poff == 1.0
        ):
            # drop the current row
            ts = ts[ts.index != ix]

    # reset the index (after removing the wrong rows)
    ts.reset_index(inplace=True)
    ts["segment_id"] = ts.index

    return ts


def sort_tss_tts_files(tss, tts, genome_size):
    """
        Sort the tss and tts files and combine them into a TS file
        # IDEA: create a function where you sort only
        and another one to merge and clean the TS file.
    """
    # combine the two file to keep the lines consistents/together
    ts = pd.merge(tss, tts, on='TUindex', how='left')

    # Added the Coding Sequence ID which will be used to predict
    # the possible transcripts that we can have
    # A transcript can contain one or more Coding Sequence
    ts["segment_id"] = ts.index

    # Cleaning the TS file
    ts = clean_sorted_ts(ts)

    # add the TSS id which gives us the ability to know whether two or more transcripts
    # are sharing the same TSS or not, this will help us if two or more picked transcripts
    # with the same TSS are select.
    ts['TSS_id'] = ts.groupby('TSS_pos').grouper.group_info[0]

    # we check if there any gene crossing the Origin
    gene_crossing_origin, last_gene_size = is_gene_crossing_origin(ts, genome_size)
    if gene_crossing_origin:
        # if it's the case it sould be moved
        # from the first line to the last one in tss/tts files
        ts = ts.apply(np.roll, shift=-1)

    return ts, gene_crossing_origin, last_gene_size


def get_TUs(ts):
    """
    Get the Transcription Units

    returns: TUs, a dictionary which have the following structure:
            {
                TU_id1: [row1, row2],
                TU_id2: [row3],
                TU_id3: [row4, row5, row6]
            }
            The TU_ids are generated and the rows are extracted from the TS file.
            Each row represent a Segment
    """
    # start grouping TUs
    TU_id = 0
    TU_rows_list = []
    # TUs = {0: [row1], 1: [row2, row3],...}
    TUs = {}

    for index, row in ts.iterrows():
        while True:
            TU_rows_list.append(row)
            TUs[TU_id] = TU_rows_list
            if row.TTS_proba_off == 1.0:
                TU_id += 1
                TU_rows_list = []
                break
            else:
                break
    return TUs


def get_tr_combinations(segment_ids, tr_id):
    """
    Get all possible transcripts
    """
    tr_combinations = {}
    s = tuple(segment_ids)
    for size in range(1, len(s) + 1):
        for index in range(len(s) + 1 - size):
            tr_combinations[tr_id] = segment_ids[index:index + size]
            tr_id += 1
            yield tr_combinations


def get_TU_tr_ids(TUs):
    """
    After get the TUs, now get the transcript id
    along with their coding segment(s) associeted with them

    returns a dictionary of the following format:
    {TU_id1: {tr_id1: [segment_id], tr_id2: [segment_id, segment_id]}

    Example:
    {0: <--- This is TU with id 0, That have 3 possible transcripts
             and each transcript is formed from one or more coding segment
        {
         0: <--- transcript 1
            [0], <--- transcript 1 contains one segment with id 0
         1: <--- transcript 2
            [1], <--- transcript 2 contains one segment with id 1
         2: <--- transcript 3
            [0, 1] <--- transcript 3 contains two segments with ids 0 and 1
        }
    }, ...
    """
    tr_id = 0
    # the possible combinations of transcripts that we can have in each TU
    TU_tr_combi = {}
    # Now, for each TU, get all combinations of possible transcripts
    for TU_id, TU_info in TUs.items():
        # get the Coding Sequences of the this TU
        segment_ids = [segment.segment_id for segment in TU_info]
        # guess all possible combinations for each TU
        tr_combinations = list(get_tr_combinations(segment_ids, tr_id))[0]
        TU_tr_combi[TU_id] = tr_combinations
        tr_id += len(TU_tr_combi[TU_id])
    return TU_tr_combi


# using TU_tr_ids we can know crate the data frame and save it in a file ;)
def get_tr_info(TU_tr_ids, ts):
    """
    Get all information of each possible transcript
    """
    columns = [
        'TUindex', 'tr_id', 'tr_TSS_pos', 'tr_TSS_strength',
        'tr_TUorient', 'tr_TTS_pos', 'tr_rate',
        'tr_segment_count', 'TSS_id'
    ]
    # create the empty dataframe
    tr_info = pd.DataFrame(columns=columns)
    # one_tr_info: will contain information about one transcript
    one_tr_info = {}

    # for each TU
    for TU_id, tr_ids_info in TU_tr_ids.items():
        # for each transcript
        for tr_id, segment_ids in tr_ids_info.items():
            # for each segment
            for segment_id in segment_ids:
                # get the row of the current segment_id
                current_row = ts.loc[ts['segment_id'] == segment_id]
                # using this (current_row) info we can build the tr_info dictionnary
                one_tr_info = {
                    'TUindex': TU_id,
                    'tr_id': tr_id,
                    'tr_TSS_pos': current_row["TSS_pos"].values[0],
                    'tr_TSS_strength': current_row["TSS_strength"].values[0],
                    'tr_TUorient': current_row["TUorient_x"].values[0],
                    'tr_TTS_pos': current_row["TTS_pos"].values[0],
                    'tr_rate': current_row["TTS_proba_off"].values[0],
                    'TSS_id': current_row["TSS_id"].values[0],
                    # tr_segment_count: how many segments are in this Transcript
                    'tr_segment_count': len(segment_ids)
                }
                # we fill the dataframe
                tr_info.loc[tr_id] = pd.Series(one_tr_info)

    # SOME CLEANING
    # This line removes duplicates, it is very important to get
    # rid of TSS<-->TSS and TTS<-->TTS regions that are considered as segment
    tr_info.drop_duplicates(['tr_TSS_pos', 'tr_TTS_pos'], inplace=True, keep='last')
    # reset the index (after removing the duplicated lines)
    tr_info.reset_index(inplace=True)
    # assign the index to tr_id
    tr_info['tr_id'] = tr_info.index.values
    # remove the index column
    tr_info.drop('index', axis=1, inplace=True)
    return tr_info


def calc_proba_off(tr_info):
    """
    Calculate the Poff probabilities
    """
    # the cursor started from the first TU
    TU_cursor = 0
    proba_rest = 1
    # loop through the list of TU IDs
    for r, TU in enumerate(tr_info["TUindex"].tolist()):
        # if we are not in the same TU (if we moved to the next one)
        if TU_cursor != TU:
            # inrement the cursor
            TU_cursor += 1
            # set the probability to 1
            proba_rest = 1

        if proba_rest > 0:
            # get the Kon and Proba_off of the current line
            current_kon = tr_info.iloc[r]['tr_TSS_strength']
            current_proba_off = tr_info.iloc[r]['tr_rate']
            # calculate the proba_off of the current transcript
            tr_proba_off = current_kon * (current_proba_off * proba_rest)
            # And update it value in the current row 'r'
            tr_info.at[r, 'tr_rate'] = tr_proba_off
            proba_rest = (1 - current_proba_off) * proba_rest
    tr_info.to_csv("transcripts_info.csv", sep="\t", index=False)
    return tr_info


def f_prob_init_rate(init_rate, sum_init_rate, DELTA_T):
    return (1 - np.exp(-sum_init_rate * DELTA_T)) * (init_rate / sum_init_rate)


def f_prob_unhooked_rate(sum_Kon, DELTA_T, RNAPs_unhooked_nbr):
    return np.exp(-sum_Kon * DELTA_T) / RNAPs_unhooked_nbr


# Get the transciption unit with the list of tts belonging to TU.
def get_TU_tts(tss, tts):
    TU_tts = col.defaultdict(list)
    for index, TUindex in enumerate(tss['TUindex'].values):
        TU_tts[TUindex].append(tts['TTS_pos'][index])
    return TU_tts


def calc_sigma(Barr_sigma, GYRASE_CONC, k_GYRASE, x0_GYRASE, GYRASE_CTE, TOPO_CONC, k_TOPO, x0_TOPO, TOPO_CTE, DELTA_T):
    d_sigma = (-GYRASE_CONC * (1 / (1 + np.exp(-k_GYRASE * (Barr_sigma - x0_GYRASE)))) * GYRASE_CTE + TOPO_CONC * 1 / (
                1 + np.exp(k_TOPO * (Barr_sigma - x0_TOPO))) * TOPO_CTE) * DELTA_T
    Barr_sigma += d_sigma

    return Barr_sigma


###################### Saving files #######################

def save_files(output_path,
               Barr_pos, Barr_type, Dom_size, Barr_ts_remain, Barr_sigma,
               tr_info, tr_nbr, tr_times, save_RNAPs_info, save_tr_info,
               save_Dom_sigma, save_Barr_pos, save_mean_sig_wholeGenome, save_Dom_size,
               DELTA_X, RNAPs_genSC,
               RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id, RNAPs_hooked_id,
               RNAPs_strand, ts_beg, ts_remain, save_nbr_RNAPs_hooked,
               init_rate, Kon, RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC):
    # make sure that the output direcory exists, and create one if it doesn't
    os.makedirs("%s/resume_sim" % output_path, exist_ok=True)
    os.makedirs("%s/all_res" % output_path, exist_ok=True)

    # add the tr_nbr to tr_info dataframe
    tr_info['number of transcripts generated'] = tr_nbr

    # remove 'tr_segment_count' and 'TSS_id' columns before saving
    tr_info = tr_info.drop(['tr_segment_count', 'TSS_id'], axis=1)

    # tr_info
    tr_info.to_csv("%s/all_tr_info.csv" % output_path, sep='\t', index=False)

    # convert tr_times dict to pandas serie
    tr_times = pd.DataFrame.from_dict(tr_times, orient='index')
    # save the tr_times to csv file
    tr_times.to_csv("%s/save_tr_times.csv" % output_path, sep='\t', index=True, header=False)
    # save the number of RNAPs hooked
    np.savez("%s/save_nbr_RNAPs_hooked.npz" % output_path, nbr_RNAPs_hooked=save_nbr_RNAPs_hooked)

    # Save last info
    np.savez("%s/resume_sim/resume_sim_RNAPs.npz" % output_path, RNAPs_tr=RNAPs_tr,
             RNAPs_pos=RNAPs_pos,
             RNAPs_unhooked_id=RNAPs_unhooked_id,
             RNAPs_strand=RNAPs_strand,
             ts_beg=ts_beg,
             ts_remain=ts_remain,
             RNAPs_hooked_id=RNAPs_hooked_id)

    np.savez("%s/resume_sim/resume_sim_tr.npz" % output_path, tr_nbr=tr_nbr,
             init_rate=init_rate)

    np.savez("%s/resume_sim/resume_sim_Barr.npz" % output_path, Barr_pos=Barr_pos,
             Barr_type=Barr_type,
             Dom_size=Dom_size,
             Barr_ts_remain=Barr_ts_remain,
             Barr_sigma=Barr_sigma)

    # Save all info
    np.savez("%s/all_res/save_RNAPs_info" % output_path, RNAPs_info=save_RNAPs_info)
    np.savez("%s/all_res/save_tr_info" % output_path, tr_info=save_tr_info)
    np.savez("%s/all_res/save_sigma_info" % output_path, dom_sigma_info=save_Dom_sigma, save_Barr_pos=save_Barr_pos,
             mean_sig_wholeGenome=save_mean_sig_wholeGenome, Dom_size=save_Dom_size)


###########################################################
#         Transcription Process (Simulation)              #
###########################################################

def start_transcribing(INI_file, output_path=None, resume=False):
    """
    Function to start/resume the transcription simulation process

    Args:
        INI_file (str): The path to the parameter file.
        output_path (str, optional): The path in which all the generated files will be saved.
        resume (bool, optional): Whether this function is used to resume an already done simulation or not.

    Returns:
        GFF_file : Relative path to the GFF file
        TSS_file : Relative path to the TSS file
        TTS_file : Relative path to the TTS file
        SIM_TIME : Simulation time
        RNAPS_NB : Numbers of RNAPol
        tr_nbr : Number of transcripts
        tr_times : Transcription time
        init_rate : Initiation rate
        RNAPs_tr : Contains the id of the picked transcript
        RNAPs_pos : The position of RNAPols
        RNAPs_unhooked_id : id of Unhoocked RNAPols
        save_RNAPs_info : Contains the RNAPs_tr and RNAPs_pos
        save_tr_info : Contains the tr_nbr and init_rate
        save_Dom_sigma : Contains the SC density (sigma) value of each domaine
        save_Barr_pos : Contains the barriers positions (whether Barr_fix or RNAPol)
        cov_bp : Coverage
        tr_end : The end (position) of transcripts
    """

    ###########################################################
    #                 initiation of variables                 #
    ###########################################################

    ####################### Params info ###################

    config = read_config_file(INI_file)

    # get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF')
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    Prot_file = config.get('INPUTS', 'BARR_FIX')

    # get values from the config file
    m = config.getfloat('PROMOTER', 'm')
    sigma_t = config.getfloat('PROMOTER', 'sigma_t')
    epsilon = config.getfloat('PROMOTER', 'epsilon')

    DELTA_X = config.getfloat('GLOBAL', 'DELTA_X')
    DELTA_T = config.getfloat('GLOBAL', 'DELTA_T')

    RNAPs_genSC = config.getfloat('SIMULATION', 'RNAPs_genSC')
    SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0')
    RNAPS_NB = config.getint('SIMULATION', 'RNAPS_NB')
    SIM_TIME = config.getfloat('SIMULATION', 'SIM_TIME')
    OUTPUT_STEP = config.getfloat('SIMULATION', 'OUTPUT_STEP')
    GYRASE_CONC = config.getfloat('SIMULATION', 'GYRASE_CONC')
    TOPO_CONC = config.getfloat('SIMULATION', 'TOPO_CONC')

    TOPO_CTE = config.getfloat('TOPOISOMERASES', 'TOPO_CTE')
    GYRASE_CTE = config.getfloat('TOPOISOMERASES', 'GYRASE_CTE')
    k_GYRASE = config.getfloat('TOPOISOMERASES', 'k_GYRASE')
    x0_GYRASE = config.getfloat('TOPOISOMERASES', 'x0_GYRASE')
    k_TOPO = config.getfloat('TOPOISOMERASES', 'k_TOPO')
    x0_TOPO = config.getfloat('TOPOISOMERASES', 'x0_TOPO')
    # Calculate SIGMA_0 based on Topoisomerases concentration.
    # SIGMA_0 = 0 #((-np.log(((GYRASE_CONC*GYRASE_CTE)/TOPO_CONC*TOPO_CTE)-1))/k)+x_0

    # path to the input files (remove the "params.ini" from the path)
    pth = INI_file.rpartition("/")[0]  # + "/"
    # or you can use : pth = os.path.split(INI_file)[0] + "/"
    # if pth=="/":
    # pth="./"
    if pth == "":
        pth = "."
    if pth[-1] != "/":
        pth += "/"

    gff_df_raw = load_gff(pth + GFF_file)

    # get the genome size
    genome_size = get_genome_size(gff_df_raw)

    tss = load_tab_file(pth + TSS_file).sort_values(by=["TSS_pos"])
    tts = load_tab_file(pth + TTS_file).sort_values(by=["TTS_pos"])
    # we will load the prot_file later

    # here we sort the TSSs and TTSs info before using them
    # and we have to take into account cases where gene is crossing the Origin
    ts, gene_crossing_origin, last_gene_size = sort_tss_tts_files(tss, tts, genome_size)

    TUs = get_TUs(ts)
    TU_tr_ids = get_TU_tr_ids(TUs)
    tr_info = get_tr_info(TU_tr_ids, ts)
    tr_info = calc_proba_off(tr_info)

    # get all possible transcripts and their info
    tr_id = tr_info["tr_id"].values
    tr_TU = tr_info["TUindex"].values
    tr_strand = str2num(tr_info["tr_TUorient"].values)
    tr_start = tr_info["tr_TSS_pos"].values
    tr_end = tr_info["tr_TTS_pos"].values
    tss_strength = tr_info["tr_TSS_strength"].values
    tr_size = abs(tr_end - tr_start)
    tr_rate = tr_info["tr_rate"].values
    tss_id = tr_info["TSS_id"].values
    ts_beg_all_trs = np.zeros(len(tr_id), dtype=int)
    ts_remain_all = tr_size

    # The RNAPs id
    RNAPs_id = np.full(RNAPS_NB, range(0, RNAPS_NB), dtype=int)
    # RNAPs_last_pos
    RNAPs_last_pos = np.full(RNAPS_NB, np.nan)

    # tr_def can e replaced by tr_info
    tr_def = pd.DataFrame(
        data={
            "tr_id": tr_id,
            "tr_TU": tr_TU,
            "tr_strand": tr_strand,
            "tr_start": tr_start,
            "tr_end": tr_end,
            "tr_size": tr_size,
            "tss_strength": tss_strength,
            "tr_rate": tr_rate
        },
        columns=["tr_id", "tr_TU", "tr_strand", "tr_start", "tr_end", "tr_size", "tss_strength", "tr_rate"])

    print("====== [BEGIN] General Information ======")
    print("genome_size ---= ", genome_size)
    print("RNAPs_id ---= ", RNAPs_id)
    print("------ Transcripts Info ------")
    print(tr_def)
    print("======= [END] General Information =======")
    #input("Press Enter to continue...")

    # Devide by DELTA_X
    genome = int(genome_size / DELTA_X)
    tr_start = (tr_start / DELTA_X).astype(int)
    tr_end = (tr_end / DELTA_X).astype(int)
    tr_size = (tr_size / DELTA_X).astype(int)
    ts_remain_all = np.around(tr_size)

    if not resume:
        # The position of RNAPs
        RNAPs_pos = np.full(RNAPS_NB, np.nan)

        # The number of times transcripts has been transcribed
        tr_nbr = np.zeros(len(tr_id), dtype=int)

        # try to read the prot_file which contains the fixed barriers positions
        try:
            prot = load_tab_file(pth + Prot_file)
            Barr_fix = (prot['prot_pos'].values / DELTA_X).astype(int)

            # just for the echo we can assign it directely
            Barr_pos = np.copy(Barr_fix)

            # update in case where no fixed barriers !!!
            # abs in case we have Barr_pos[i+1]>Barr_pos[i] e.g: [64 57]
            Dom_size = np.abs(np.ediff1d(Barr_pos))  # , dtype=int

            Dom_size = np.append(Dom_size, genome - (np.max(Barr_pos) - np.min(Barr_pos)))

            # Barr_type contains the barrier type
            # 0 : fixed barrrier (e.g Protein)
            # -1 : -RNAPol (direction : <--)
            # 1 : +RNAPol (direction : -->)
            Barr_type = np.full(len(Barr_fix), 0, dtype=int)
            Barr_sigma = np.full(len(Barr_fix), SIGMA_0)
            # here we need to create Barr_ts_remain
            # to track the position of each RNAPol
            # each position in Barr_ts_remain is associated with the same position in Barr_pos
            Barr_ts_remain = np.full(len(Barr_fix), np.nan)  # The Barr_ts_remain of fixed barr is NaN

        # if prot_file is empty or doesn't exist then:
        except (pd.io.common.EmptyDataError, OSError, ValueError, KeyError):
            # we'll have one Dom_size which is the whole genome
            # There is no Barr_fix
            Dom_size = np.array([genome], dtype=int)
            # create the other variables
            Barr_type = np.array([], dtype=int)
            Barr_sigma = np.array([SIGMA_0])  # on the whole genome
            Barr_pos = np.array([], dtype=int)
            Barr_ts_remain = np.array([])

        RNAPs_unhooked_id = np.copy(RNAPs_id)
        RNAPs_strand = np.full(RNAPS_NB, np.nan)
        ts_beg = np.full(RNAPS_NB, np.nan)
        ts_remain = np.full(RNAPS_NB, np.nan)
        # RNAPs_tr will contain the id of the picked transcript
        RNAPs_tr = np.full(RNAPS_NB, -1, dtype=(int))
        # get the TSSs ids
        ts_id = ts.index.values

        # in the case of RNAP_NBR = 0
        RNAPs_hooked_id = []

        if output_path is None:
            # no resume ==> set the 'output_path' to 'first_output'
            output_path = "output"
            # 'output_path' variable is used to save the output files
            # in 'save_files()' function

    else:
        # RNAPs_info contains : ['RNAPs_unhooked_id', 'RNAPs_pos', 'RNAPs_tr']
        RNAPs_info = np.load("%s/resume_sim/resume_sim_RNAPs.npz" % output_path)

        # get the RNAPs position
        RNAPs_pos = RNAPs_info['RNAPs_pos']

        # The number of times transcripts has been transcribed
        csv_path = "%s/save_tr_nbr.csv" % output_path
        tr_nbr = get_tr_nbr_csv(csv_path)

        # when we resume, we won't need info from 'prot' file because we already
        # have all what we need in 'resume_sim_Barr.npz' file ;)

        # Get info from NPZ file
        Barr_info = np.load("%s/resume_sim/resume_sim_Barr.npz" % output_path)

        Barr_pos = Barr_info['Barr_pos']
        Dom_size = Barr_info['Dom_size']

        Barr_type = Barr_info['Barr_type']
        Barr_sigma = Barr_info['Barr_sigma']

        # here we need to make an Barr_ts_remain
        # so we can track the position of each RNAPol
        # each position in Barr_ts_remain is associated with the same position in Barr_pos
        Barr_ts_remain = Barr_info['Barr_ts_remain']

        ## do the same for RNAPs_info
        # get the RNAPs_info
        RNAPs_info = np.load("%s/resume_sim/resume_sim_RNAPs.npz" % output_path)

        # get the RNAPs_hooked_id and RNAPs_pos
        RNAPs_hooked_id = RNAPs_info["RNAPs_hooked_id"]
        RNAPs_pos = RNAPs_info["RNAPs_pos"]
        # deduce the RNAPs_hooked_pos from the extracted info
        RNAPs_hooked_pos = RNAPs_pos[RNAPs_hooked_id].astype(int)

        # since we continue the simulation, we shall retrieve the RNAPs_unhooked_id from the npz file
        RNAPs_unhooked_id = RNAPs_info["RNAPs_unhooked_id"]

        RNAPs_strand = RNAPs_info["RNAPs_strand"]
        ts_beg = RNAPs_info["ts_beg"]
        ts_remain = RNAPs_info["ts_remain"]
        # RNAPs_tr contains the id of the picked transcript
        RNAPs_tr = RNAPs_info["RNAPs_tr"]
        # get the TSSs ids
        ts_id = ts.index.values

        # will do the same for RNAPs_hooked_id
        RNAPs_hooked_id = RNAPs_info["RNAPs_hooked_id"]

        if output_path is None:
            # resume ==> set the 'output_path' to 'resume_output'
            output_path = "resume_output"

    ######### Variables used to get the coverage ##########

    id_shift_fwd = list(range(1, genome))
    id_shift_fwd.append(0)
    id_shift_fwd = np.array(id_shift_fwd)
    id_shift_bwd = list(range(0, genome - 1))
    id_shift_bwd.insert(0, genome - 1)
    id_shift_bwd = np.array(id_shift_bwd)

    cov_bp = np.arange(0, genome_size, DELTA_X)
    cov_bp = np.resize(cov_bp, genome)

    # save the time when RNApoly is starting trasncribing a specific transcript
    # tr_times = col.defaultdict(list)
    tr_times = {}
    for transcript in tr_id:
        tr_times[transcript] = []

    # numpy array where all RNAPs info will be saved except the nbr_RNAPs_hooked
    save_RNAPs_info = np.full([RNAPS_NB, 2, int(SIM_TIME / (DELTA_T * OUTPUT_STEP))], np.nan)  # nbr d'ele (cols)

    # this array will contain the number of RNAPs hooked at each time step
    # the length of the array is equivalent to the simulation length
    save_nbr_RNAPs_hooked = np.full(int(SIM_TIME / (DELTA_T * OUTPUT_STEP)), np.nan)

    # the same for transcripts info
    save_tr_info = np.full([len(tr_id), 2, int(SIM_TIME / (DELTA_T * OUTPUT_STEP))], np.nan)

    # in those variables, we will save/append info in each time step to save them as --> all_res
    save_Dom_sigma = list()
    save_Dom_size = list()
    save_Barr_pos = list()
    save_mean_sig_wholeGenome = list()
    save_Barr_pos = list()

    # The main loop
    for t in range(0, int(SIM_TIME / DELTA_T)):
        # we need to know each TSS belong to which Domaine
        # if we have only one domaine (e.g no Barr_fix) then
        # all the TSS belong to the only domaine that exists (e.g TSS_pos_idx will be 0)
        TSS_pos_idx = np.searchsorted(Barr_pos, tr_start)

        # after knowing the domaine of each TSS we can get sigma
        try:
            sigma_tr_start = Barr_sigma[TSS_pos_idx - 1]
        except IndexError:
            sigma_tr_start = np.array([SIGMA_0])

        # get the initiation rates
        # calculate the initiation rate of each transcript/gene
        init_rate = f_init_rate(tss_strength, sigma_tr_start, sigma_t, epsilon, m)
        # use the calculated init_rate to get the probability
        # of RNAPol's binding to each TSS (e.g prob_init_rate)
        sum_init_rate = np.sum(init_rate)
        prob_init_rate = f_prob_init_rate(init_rate, sum_init_rate, DELTA_T)

        if np.size(RNAPs_unhooked_id) != 0:
            # get the probability of RNAPol's stay unhooked
            prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rate, DELTA_T, len(RNAPs_unhooked_id))
            # craete the numpy array
            prob_unhooked_rate = np.full(len(RNAPs_unhooked_id), prob_unhooked_rate)
            all_prob = np.concatenate([prob_init_rate, prob_unhooked_rate])
            # create the numpy array that will contains [ nTSS , Unhooked RNAPS ]
            # e.g if we have 2 TSSs and 3 unhooked RNAPols then
            # tss_and_unhooked_RNAPs = [0, 1, -1, -1, -1]
            tss_and_unhooked_RNAPs = np.concatenate([ts_id, np.full(len(RNAPs_unhooked_id), -1, dtype=int)])
            # pick up a random transcipt
            picked_tr = np.random.choice(tss_and_unhooked_RNAPs, len(RNAPs_unhooked_id), replace=False,
                                         p=all_prob)
            # This is the KEY !
            picked_tr_hooked_id = picked_tr[np.where(picked_tr != -1)[0]]
            picked_tr_unhooked_id = picked_tr[np.where(picked_tr == -1)[0]]

            # if tss_id is duplicated
            # TODO: replace/remove tss_id with tr_start?
            unique_tss_id, counts = np.unique(tss_id, return_counts=True)
            # if more than one transcript is picked
            # and there are transcripts that share the same TSS
            if np.size(picked_tr_hooked_id) > 1 and np.size(unique_tss_id) < np.size(tss_id):
                filtred_picked_tr_hooked_id = np.array([])
                # the picked tss id
                picked_tss_id = tr_info['TSS_id'][picked_tr_hooked_id].values
                # for each tss
                for tss in np.unique(picked_tss_id):
                    # get all tr_id belonging to this tss
                    all_tss_tr_id = tr_info[tr_info['TSS_id'] == tss]['tr_id'].values
                    # from all_tss_tr_id, select only the ones that are picked
                    mask_this_tss_picked_tr_hooked_id = np.isin(all_tss_tr_id, picked_tr_hooked_id)                    
                    this_tss_picked_tr_hooked_id = all_tss_tr_id[mask_this_tss_picked_tr_hooked_id]
                    # all the transcripts have equa probability of being chosen
                    equal_prob = np.repeat(1/np.size(this_tss_picked_tr_hooked_id), np.size(this_tss_picked_tr_hooked_id))
                    chosen_tr = np.random.choice(this_tss_picked_tr_hooked_id, 1, replace=False, p=equal_prob)
                    filtred_picked_tr_hooked_id = np.append(filtred_picked_tr_hooked_id, chosen_tr).astype(int)

                # IDEA: get index and use delete may be faster
                filtred_out_picked_tr_hooked_id = np.setdiff1d(picked_tr_hooked_id, filtred_picked_tr_hooked_id)
                picked_tr[np.isin(picked_tr, filtred_out_picked_tr_hooked_id)] = -1
                # reset it again!
                picked_tr_hooked_id = picked_tr[np.where(picked_tr != -1)[0]]

            new_RNAPs_hooked_id = RNAPs_unhooked_id[np.where(picked_tr != -1)[0]]

            RNAPs_tr[new_RNAPs_hooked_id] = picked_tr[picked_tr != -1]
            RNAPs_strand[new_RNAPs_hooked_id] = tr_strand[picked_tr[np.where(picked_tr != -1)]].astype(int)

            # The new position of each polymerase
            # if there is no RNAP already at this position
            RNAPs_pos[new_RNAPs_hooked_id] = tr_start[picked_tr[np.where(picked_tr != -1)]].astype(int)

            # Bug #1 Fix
            # We use new_RNAPs_hooked_id and RNAPs_pos[new_RNAPs_hooked_id] to create an ordered dictionary
            # new_hooked_RNAPs_pos_dict contains [new hooked RNAP position as a key: new RNAP hooked id as value]
            new_hooked_RNAPs_pos_dict = dict(zip(RNAPs_pos[new_RNAPs_hooked_id], new_RNAPs_hooked_id))
            # We sort the dict by Keys (e.g position)
            # At the end we'll get the sorted positions with their corresponding ids
            new_hooked_RNAPs_pos_ordered_dict = col.OrderedDict(sorted(new_hooked_RNAPs_pos_dict.items()))
            # Here idx are sorted based on the positions
            new_hooked_RNAPs_idx_sorted = np.array([idx for idx in new_hooked_RNAPs_pos_ordered_dict.values()]).astype(int)
            new_hooked_RNAPs_pos_sorted = np.array([pos for pos in new_hooked_RNAPs_pos_ordered_dict.keys()]).astype(int)

            # take the positions and use them to get the index in which we will insert
            # the new recruited RNAPs in Barr_pos array
            Barr_pos_RNAPs_idx = np.searchsorted(Barr_pos, new_hooked_RNAPs_pos_sorted)
            #!print("new_hooked_RNAPs_pos_sorted ===> ", new_hooked_RNAPs_pos_sorted)
            # After we got everything we need, and the new RNAPol is ready be hooked
            # we should first check if the hooked position is empty and there is no
            # other RNAPol passing by/in the same position
            # if it's the case
            # we abord the insertion of this specific RNAPol
            next_barr_pos = np.copy(Barr_pos)
            # IDEA: Can't we deduplicate this ?
            next_barr_pos[np.where(Barr_type == -1)] -= 1
            next_barr_pos[np.where(Barr_type == 1)] += 1
            is_already_filled_pos = np.where(np.isin(next_barr_pos, new_hooked_RNAPs_pos_sorted))
            already_filled_pos = next_barr_pos[is_already_filled_pos]
            # next get the idx position of the cancelled insert new_hooked_RNAPs_pos
            already_filled_pos_idx = np.where(new_hooked_RNAPs_pos_sorted == already_filled_pos)[0]
            # after getting the index of the cancelled insertion position we remove them from Barr_pos_RNAPs_idx
            Barr_pos_RNAPs_idx = np.delete(Barr_pos_RNAPs_idx, already_filled_pos_idx)
            new_hooked_RNAPs_pos_sorted = np.delete(new_hooked_RNAPs_pos_sorted, already_filled_pos_idx)
            new_hooked_RNAPs_idx_sorted = np.delete(new_hooked_RNAPs_idx_sorted, already_filled_pos_idx)
            # get the value picked_tr_hooked_id_to_cancel
            picked_tr_hooked_id_to_cancel = picked_tr_hooked_id[already_filled_pos_idx]
            # remove them from picked_tr_hooked_id
            picked_tr_hooked_id = np.delete(picked_tr_hooked_id, already_filled_pos_idx)
            # then set them to -1 in picked_tr again
            picked_tr[np.where(picked_tr == picked_tr_hooked_id_to_cancel)[0]] = -1

            #!print("Barr_pos_RNAPs_idx =========> ", Barr_pos_RNAPs_idx)
            #!print("RNAPs_strand =========> ", RNAPs_strand)
            #!print("new_hooked_RNAPs_idx_sorted =========> ", new_hooked_RNAPs_idx_sorted)

            # if we have one or no barrier on the genome
            if Barr_pos.size <= 1:
                # if there is no new RNAPol hooked
                if new_hooked_RNAPs_idx_sorted.size == 0:
                    Dom_size = np.array([genome])
                    Barr_sigma = np.array([SIGMA_0])
                # if there is ONE or MORE new RNAPol hooked
                elif new_hooked_RNAPs_idx_sorted.size >= 1:
                    Barr_pos = np.insert(Barr_pos, Barr_pos_RNAPs_idx, new_hooked_RNAPs_pos_sorted)
                    Barr_type = np.insert(Barr_type, Barr_pos_RNAPs_idx, RNAPs_strand[new_hooked_RNAPs_idx_sorted])
                    Dom_size = np.abs(np.ediff1d(Barr_pos))
                    Dom_size = np.append(Dom_size, genome - (np.max(Barr_pos) - np.min(Barr_pos)))
                    Barr_sigma = np.repeat(SIGMA_0, np.size(Barr_pos))

            # if we have at least two barrier
            elif Barr_pos.size >= 2:
                Barr_pos = np.insert(Barr_pos, Barr_pos_RNAPs_idx, new_hooked_RNAPs_pos_sorted)
                Barr_type = np.insert(Barr_type, Barr_pos_RNAPs_idx, RNAPs_strand[new_hooked_RNAPs_idx_sorted])
                Dom_size = np.abs(np.ediff1d(Barr_pos))
                Dom_size = np.append(Dom_size, genome - (np.max(Barr_pos) - np.min(Barr_pos)))
                Barr_sigma = np.insert(Barr_sigma, Barr_pos_RNAPs_idx, Barr_sigma[Barr_pos_RNAPs_idx - 1])

            # RNAPs_last_pos
            RNAPs_last_pos[new_hooked_RNAPs_idx_sorted] = tr_end[picked_tr_hooked_id]

            ts_beg[new_hooked_RNAPs_idx_sorted] = 0
            ts_remain[new_hooked_RNAPs_idx_sorted] = ts_remain_all[picked_tr_hooked_id]  # NOT picked_tr
            #!print("picked_tr_hooked_id -----> ", picked_tr_hooked_id)
            #!print("ts_remain -----> ", ts_remain)
            #!print("ts_remain_all -----> ", ts_remain_all)
            #!print("tr_size -----> ", tr_size)
            Barr_ts_remain = np.insert(Barr_ts_remain, Barr_pos_RNAPs_idx, ts_remain[new_hooked_RNAPs_idx_sorted])
            RNAPs_hooked_id = np.where(RNAPs_tr != -1)[0]

        ts_beg[RNAPs_hooked_id] += 1
        ts_remain[RNAPs_hooked_id] -= 1

        # save the time when RNAPol FINISHS trasncribing a specific transcript
        for x in RNAPs_tr[np.where(ts_remain == 0)]:
            tr_times[x].append(t * DELTA_T)  # + 0.5

        tr_nbr[RNAPs_tr[np.where(ts_remain == 0)]] += 1

        # look in the net : numpy where two conditions
        Barr_ts_remain[np.where(Barr_type == -1)] -= 1
        Barr_ts_remain[np.where(Barr_type == 1)] -= 1

        # Get the index of RNAPs to remove
        rm_RNAPs_idx = np.where(Barr_ts_remain == 0)[0]

        # recover sigma value of the removed position
        removed_sigma = Barr_sigma[rm_RNAPs_idx]
        removed_dom_size = Dom_size[rm_RNAPs_idx]

        # recover the old_dom_size : the size of the previous domaine before combination/merging
        old_dom_size = Dom_size[rm_RNAPs_idx - 1]
        old_sigma = Barr_sigma[rm_RNAPs_idx - 1]

        # update Dom_size
        Dom_size[rm_RNAPs_idx - 1] += removed_dom_size
        # or
        # abs in case we have Barr_pos[i+1]>Barr_pos[i] e.g: [64 57]
        # Dom_size = np.abs(np.ediff1d(Barr_pos))
        # Dom_size = np.append(Dom_size, genome-Barr_fix[-1]+Barr_fix[0])

        Barr_sigma[rm_RNAPs_idx - 1] = (old_dom_size * old_sigma + removed_dom_size * removed_sigma) / (
                    old_dom_size + removed_dom_size)

        # and reomve them
        Barr_pos = np.delete(Barr_pos, rm_RNAPs_idx)
        Barr_type = np.delete(Barr_type, rm_RNAPs_idx)
        Barr_ts_remain = np.delete(Barr_ts_remain, rm_RNAPs_idx)
        Barr_sigma = np.delete(Barr_sigma, rm_RNAPs_idx)
        Dom_size = np.delete(Dom_size, rm_RNAPs_idx)

        # update the RNAPs_tr array
        RNAPs_tr[np.where(ts_remain == 0)] = -1
        # update the RNAPs_unhooked_id based on RNAPs_tr
        RNAPs_unhooked_id = np.where(RNAPs_tr == -1)[0]

        # reset the arrays
        RNAPs_strand[RNAPs_unhooked_id] = np.nan
        RNAPs_pos[RNAPs_unhooked_id] = np.nan
        RNAPs_last_pos[RNAPs_unhooked_id] = np.nan
        ts_beg[RNAPs_unhooked_id] = np.nan
        ts_remain[RNAPs_unhooked_id] = np.nan

        Barr_pos[np.where(Barr_type == -1)] -= 1
        Barr_pos[np.where(Barr_type == 1)] += 1
        # Barr_pos=np.mod(Barr_pos,genome)

        # Update the position of polymerases still transcribing
        RNAPs_pos[np.where(RNAPs_strand == 1)] += 1
        RNAPs_pos[np.where(RNAPs_strand == -1)] -= 1

        # Update the Dom_size (+1 or -1)
        # if we have at least two barrier
        try:
            # abs in case we have Barr_pos[i+1]>Barr_pos[i] e.g: [64 57]
            Dom_size = np.abs(np.ediff1d(Barr_pos))
            Dom_size = np.append(Dom_size, genome - (np.max(Barr_pos) - np.min(Barr_pos)))
        # in case we have one or zero barrier
        except (IndexError, ValueError):
            Dom_size = np.array([genome])
            #Barr_sigma = np.repeat(SIGMA_0, np.size(new_hooked_RNAPs_idx_sorted))
            Barr_sigma = np.array([SIGMA_0])

        # UPDATE SIGMA
        # R_plus_pos : the ids of RNA pol in the + strand
        R_plus_pos = np.where(Barr_type == 1)[0].astype(int)
        # R_minus_pos : the ids of RNA pol in the - strand
        R_minus_pos = np.where(Barr_type == -1)[0].astype(int)
        #### Extract all types of domaines (Those are ids of domaines)
        # Barr_type_ahead to make the extraction circular ;)
        Barr_type_ahead = np.roll(Barr_type, -1)
        # __|__________O+____
        Barr_Dom_RPlus = np.where((Barr_type == 0) & (Barr_type_ahead == 1))
        # __|__________O-____
        Barr_Dom_RMinus = np.where((Barr_type == 0) & (Barr_type_ahead == -1))
        # __|__________|_____
        Barr_Dom_Barr = np.where((Barr_type == 0) & (Barr_type_ahead == 0))
        # ___O+_________O+___
        RPlus_Dom_RPlus = np.where((Barr_type == 1) & (Barr_type_ahead == 1))
        # ___O-_________O-___
        RMinus_Dom_RMinus = np.where((Barr_type == -1) & (Barr_type_ahead == -1))
        # ___O+_________O-___
        RPlus_Dom_RMinus = np.where((Barr_type == 1) & (Barr_type_ahead == -1))
        # ___O-_________O+___
        RMinus_Dom_RPlus = np.where((Barr_type == -1) & (Barr_type_ahead == +1))
        # ___O-_________|____
        RMinus_Dom_Barr = np.where((Barr_type == -1) & (Barr_type_ahead == 0))
        # ___O+_________|____
        RPlus_Dom_Barr = np.where((Barr_type_ahead == 0) & (Barr_type == +1))

        """
        This sim_info is just for the print (can be removed)
        """
        if Barr_sigma.size == Barr_type.size:
            sim_info = pd.DataFrame(
                data={
                    "Barr_sigma": Barr_sigma,
                    "Barr_type": Barr_type,
                    "Barr_pos": Barr_pos,
                    "Dom_size": Dom_size,
                    "Barr_ts_remain": Barr_ts_remain},
                columns=["Barr_sigma", "Barr_type", "Barr_pos", "Dom_size", "Barr_ts_remain"])
            #!print(sim_info)

        #!print("======  ======  ======  ======  ======")
        #!print("======  ======  ======  ======  ======")

        #### And then correct the value of Sigma in each case (before/after)
        corr_sig_Barr_Dom_RPlus = (Dom_size[Barr_Dom_RPlus] - 1) / (Dom_size[Barr_Dom_RPlus])  # Sigma decrease x1
        corr_sig_Barr_Dom_RMinus = (Dom_size[Barr_Dom_RMinus] + 1) / (Dom_size[Barr_Dom_RMinus])  # Sigma increase x1
        corr_sig_Barr_Dom_Barr = (Dom_size[Barr_Dom_Barr]) / (Dom_size[Barr_Dom_Barr])  # Sigma FIX
        corr_sig_RPlus_Dom_RPlus = (Dom_size[RPlus_Dom_RPlus]) / (Dom_size[RPlus_Dom_RPlus])  # Sigma FIX
        corr_sig_RMinus_Dom_RMinus = (Dom_size[RMinus_Dom_RMinus]) / (Dom_size[RMinus_Dom_RMinus])  # Sigma FIX
        corr_sig_RPlus_Dom_RMinus = (Dom_size[RPlus_Dom_RMinus] + 2) / (Dom_size[RPlus_Dom_RMinus])  # Sigma increase x2
        corr_sig_RMinus_Dom_RPlus = (Dom_size[RMinus_Dom_RPlus] - 2) / (Dom_size[RMinus_Dom_RPlus])  # Sigma decrease x2
        corr_sig_RMinus_Dom_Barr = (Dom_size[RMinus_Dom_Barr] - 1) / (Dom_size[RMinus_Dom_Barr])  # Sigma decrease x1
        corr_sig_RPlus_Dom_Barr = (Dom_size[RPlus_Dom_Barr] + 1) / (Dom_size[RPlus_Dom_Barr])  # Sigma increase x1

        ### Multiply Sigma *= Corr (Each sigma value correspond to an specific domaine)
        Barr_sigma[Barr_Dom_RPlus] *= corr_sig_Barr_Dom_RPlus
        Barr_sigma[Barr_Dom_RMinus] *= corr_sig_Barr_Dom_RMinus
        Barr_sigma[Barr_Dom_Barr] *= corr_sig_Barr_Dom_Barr
        Barr_sigma[RPlus_Dom_RPlus] *= corr_sig_RPlus_Dom_RPlus
        Barr_sigma[RMinus_Dom_RMinus] *= corr_sig_RMinus_Dom_RMinus
        Barr_sigma[RPlus_Dom_RMinus] *= corr_sig_RPlus_Dom_RMinus
        Barr_sigma[RMinus_Dom_RPlus] *= corr_sig_RMinus_Dom_RPlus
        Barr_sigma[RMinus_Dom_Barr] *= corr_sig_RMinus_Dom_Barr
        Barr_sigma[RPlus_Dom_Barr] *= corr_sig_RPlus_Dom_Barr

        ### calculate the Supercoiling generated in each domaine
        # RNAPs_genSC_all : contains an array of RNAPs_genSC that should be added or substracted from each domaine
        RNAPs_genSC_all = RNAPs_genSC / Dom_size
        # update the value of sigma
        Barr_sigma[Barr_Dom_RPlus] -= RNAPs_genSC_all[Barr_Dom_RPlus]
        Barr_sigma[Barr_Dom_RMinus] += RNAPs_genSC_all[Barr_Dom_RMinus]
        Barr_sigma[RPlus_Dom_RMinus] += 2 * RNAPs_genSC_all[RPlus_Dom_RMinus]
        Barr_sigma[RMinus_Dom_RPlus] -= 2 * RNAPs_genSC_all[RMinus_Dom_RPlus]
        Barr_sigma[RMinus_Dom_Barr] -= RNAPs_genSC_all[RMinus_Dom_Barr]
        Barr_sigma[RPlus_Dom_Barr] += RNAPs_genSC_all[RPlus_Dom_Barr]
        # We shall consider the case in which we'll have one RNAPol
        # transcribing from right or the left and without any existing barrier on the other side
        # i : you can relace 'if' with None_Dom_RPlus = np.where((len(Barr_type)==1) & (Barr_type_ahead==1))
        if len(Barr_type) == 1:
            # ____________O+_____
            None_Dom_RPlus = np.where(Barr_type_ahead == 1)
            # ____________O-_____
            None_Dom_RMinus = np.where(Barr_type_ahead == -1)
            # __O+_______________
            RPlus_Dom_None = np.where(Barr_type_ahead == 1)
            # __O-_______________
            RMinus_Dom_None = np.where(Barr_type_ahead == -1)
            # update the value of sigma (one Barr case)
            Barr_sigma[None_Dom_RPlus] -= RNAPs_genSC_all[None_Dom_RPlus]
            Barr_sigma[None_Dom_RMinus] += RNAPs_genSC_all[None_Dom_RMinus]
            Barr_sigma[RPlus_Dom_None] += RNAPs_genSC_all[RPlus_Dom_None]
            Barr_sigma[None_Dom_RMinus] -= RNAPs_genSC_all[None_Dom_RMinus]

        # Now calc_sigma
        Barr_sigma = calc_sigma(Barr_sigma, GYRASE_CONC, k_GYRASE, x0_GYRASE, GYRASE_CTE, TOPO_CONC, k_TOPO, x0_TOPO,
                                TOPO_CTE, DELTA_T)

        try:
            mean_sig_wholeGenome = np.sum(Barr_sigma * Dom_size) / genome
        except ValueError:
            # in case we have one RNAPol transcribing (No Barr_fix)
            # ________O______________
            mean_sig_wholeGenome = (Barr_sigma[0] + Barr_sigma[1]) / 2

        # Update the initiation rate
        init_rate = f_init_rate(tss_strength, sigma_tr_start, sigma_t, epsilon, m)

        if t % OUTPUT_STEP == 0:
            tt = int(t // OUTPUT_STEP)
            # save all informations to npz file
            # RNAPs_info
            save_RNAPs_info[:, 0, tt] = RNAPs_tr
            save_RNAPs_info[:, 1, tt] = RNAPs_pos
            # save the number of hooked RNAPs
            save_nbr_RNAPs_hooked[tt] = np.size(RNAPs_hooked_id)
            # tr_info
            save_tr_info[:, 0, tt] = tr_nbr
            save_tr_info[:, 1, tt] = init_rate

        save_Dom_sigma.append(Barr_sigma)
        save_Dom_size.append(Dom_size)
        save_Barr_pos.append(Barr_pos)
        save_mean_sig_wholeGenome.append(mean_sig_wholeGenome)

    save_Dom_sigma = np.array(save_Dom_sigma)
    save_Dom_size = np.array(save_Dom_size)
    save_Barr_pos = np.array(save_Barr_pos)
    save_mean_sig_wholeGenome = np.array(save_mean_sig_wholeGenome)

    # make sure that the output direcory exists, and create one if it doesn't
    os.makedirs(output_path, exist_ok=True)

    try:
        # Copy the params to the output folder
        copy(INI_file, output_path)
    except Exception as e:
        print("Input file was not copied")
        sys.exit(1)

    save_files(output_path, Barr_pos, Barr_type, Dom_size, Barr_ts_remain, Barr_sigma, tr_info, tr_nbr, tr_times,
               save_RNAPs_info, save_tr_info, save_Dom_sigma, save_Barr_pos, save_mean_sig_wholeGenome, save_Dom_size,
               DELTA_X, RNAPs_genSC, RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id, RNAPs_hooked_id, RNAPs_strand, ts_beg,
               ts_remain, save_nbr_RNAPs_hooked, init_rate, tr_rate, RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC)

    print("Simulation completed successfully !! \nNumber of transcripts :")
    for i, v in enumerate(tr_nbr):
        print("Transcript ID {} : {}".format(i, v))

    return (GFF_file, TSS_file, TTS_file,
            SIM_TIME, RNAPS_NB,
            tr_nbr, tr_times, init_rate,
            RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id,
            save_RNAPs_info, save_tr_info, save_Dom_sigma, save_Barr_pos,
            cov_bp, tr_end)
