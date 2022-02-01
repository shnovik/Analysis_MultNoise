# Description: The purpose of this file is to be able to take in start/stop times data for Up states, Down states, and SWS periods in one recording, and calculate statistics

# python modules
import numpy as np
import math
from scipy import fft
from scipy.optimize import curve_fit
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib import cycler

import sys
import csv

# Statistics calculations:

def mean_cv_skew(up_dur,down_dur):
    # This function taken from Ann Christin's Masterthesis code
    """Calculate the mean, CV and skew of the durations"""
    UPdur_array = np.array(up_dur)
    DOWNdur_array = np.array(down_dur)

    mean_up = np.mean(UPdur_array)
    std_up = np.std(UPdur_array)
    mean_third_up = np.mean(UPdur_array**3)

    cv_up = std_up/mean_up
    skew_up = (mean_third_up-3*mean_up*(std_up**2)-mean_up**3)/(std_up**3)

    mean_down = np.mean(DOWNdur_array)
    std_down = np.std(DOWNdur_array)
    mean_third_down = np.mean(DOWNdur_array**3)

    cv_down = std_down/mean_down
    skew_down = (mean_third_down-3*mean_down*(std_down**2)-mean_down**3)/(std_down**3)

    return mean_up,cv_up,mean_down,cv_down,skew_up,skew_down


def SCC_many_swsperiods(U, D, lag):
    """
    For lag = 0 or 1, returns a numpy array of length 1, containing either the U->D SCC or the D->U SCC

    Inputs:
        U: A list of lists, where each list contains the Up durations in one SWS period from one recording
        D: A list of lists, where earch list contains the Down durations in one SWS period from one recording
        lag: lag = 0 for calculationg D->U SCC, lag = 1 for calculating U->D SCC

    Outputs: The desired SCC
    """

    if lag < 0:
        X = D
        Y = U
        lag= -lag
    else:
        X = U
        Y = D
    
    # Compute the means and standard deviations
    used_X = []
    used_Y = []
    for swsIndex in range(len(X)):
        L = len(X[swsIndex])
        part_used_X = X[swsIndex][:L-lag]
        part_used_Y = Y[swsIndex][lag:]
        used_X.extend(part_used_X)
        used_Y.extend(part_used_Y)
    used_X_array = np.array(used_X, dtype = object)
    used_Y_array = np.array(used_Y, dtype = object)
    
    m1=np.mean(used_X_array)
    m2=np.mean(used_Y_array)
    s1=np.std(used_X_array)
    s2=np.std(used_Y_array)


    # Calculate covariance using all truly consecutive pairs (pairs in the same sws period)
    sws_covariance = []
    for swsIndex in range(len(X)):
        swsX = np.array(X[swsIndex])
        swsY = np.array(Y[swsIndex])
        L = len(swsX)
        sum_terms = (swsX[:L-lag] - m1) * (swsY[lag:] - m2)
        sws_covariance.extend(sum_terms)
    C = np.mean(sws_covariance)
    return C/(s1*s2)



def SCC(U,D,lag):
    # This function taken from Ann Christin's Masterthesis code
    """
    returns vector of length lag. lag is a number so we return the correlation (float) for this lag
    of serial correlation coefficients of stationary sequence X
    according to definition in Jercog et al
    We assume for this function that we start in Down. D_i,U_i,D_i+1...
    lag <= 0 refers to D before U
    In particular lag = 0 refers to the immediately previous DOWN.
    lag = 1 refers to the immediately consecutive Down.
    lag > 0 refers to DOWN after UP
    """

    L=len(U)
    if lag<0:
        X=D
        Y=U
        lag= -lag
    else:
        X=U
        Y=D

    m1=np.mean(X[:L-lag])
    m2=np.mean(Y[lag:])
    #print("mean Down/Up",m1)
    #print("mean Down/Up",m2)
    s1=np.std(X[:L-lag])
    s2=np.std(Y[lag:])
    #print("standard dev s1",s1)
    #print("standard dev s2",s2)

    #first = np.subtract(X[:L-lag], m1)
    #second = np.subtract(Y[lag:], m2)
    #sumTerms = np.multiply(first, second[0:len(first)])
    #C = np.mean(sumTerms)

    C=np.mean(X[:L-lag]*Y[lag:])-m1*m2

    return C/(s1*s2)


def ratio_simtime_up_down(up_dur,down_dur):
    # This function taken from Ann Christin's Masterthesis code
    """Compute ratio of simulation time spent in Up or Down state
    For the denominator we consider the up + down durations and not the whole simulation time because parts of the simulation time are not included."""

    time_in_up = np.sum(up_dur)
    time_in_down = np.sum(down_dur)
    total_time = time_in_up + time_in_down

    perc_up = np.round(time_in_up/total_time*100,2)
    perc_down = np.round(time_in_down/total_time*100,2)

    return perc_up,perc_down


def calculateStats(up_file, down_file, sws_file):
    """Given the Up state start/stop times, the Down state start/stop times, and the overall SWS period start/stop times for one recording, calculate 8 statistics
    
    Inputs:
        up_file   -- file name of a csv file containing start and stop times for many up states
        down_file -- file name of a csv file containing start and stop times for many down states
        sws_file  -- file name of a csv file containing start and stop times of slow wave sleep episodes

    Outputs: (all are floats)
        mean_up, cv_up, mean_down, cv_down, ud, du, perc_up, perc_down
    
    """
    # Create arrays of up and down durations
    up_times_list = list(csv.reader(open(up_file)))
    down_times_list = list(csv.reader(open(down_file)))
    sws_times_list = list(csv.reader(open(sws_file)))
    up_times = np.array(up_times_list)
    down_times = np.array(down_times_list)
    sws_times = np.array(sws_times_list)

    # Initialize lists for storing statistics
    Uscc = []
    Dscc = []
    Umcv = []
    Dmcv = []

    # Collect the durations in lists
    for swEpisode in sws_times[1:]:
        sws_start = float(swEpisode[0])
        sws_stop = float(swEpisode[1])

        up_dur = []
        down_dur = []
        for timePair in up_times[1:]:
            state_start = float(timePair[0])
            state_stop = float(timePair[1])
            if state_start >= sws_start and state_stop <= sws_stop:
                dur = state_stop - state_start
                up_dur.append(dur)
        for timePair in down_times[1:]:
            state_start = float(timePair[0])
            state_stop = float(timePair[1])
            if state_start >= sws_start and state_stop <= sws_stop:
                dur = state_stop - state_start
                down_dur.append(dur)
        
        
        if len(up_dur) == len(down_dur) - 1: #this is typical, since all SWS periods start with a D, and many also end with a D
            down_dur = down_dur[:-1]

        if len(up_dur) == len(down_dur): #skip that one SWS period that somehow has 2 more Downs than Ups
            Uscc.append(up_dur)
            Dscc.append(down_dur)
        

        Umcv.extend(up_dur)
        Dmcv.extend(down_dur)   

    # Calculate statistics
    mean_up, cv_up, mean_down, cv_down, skew_up, skew_down = mean_cv_skew(Umcv, Dmcv)
    perc_up, perc_down = ratio_simtime_up_down(Umcv, Dmcv)
    ud = SCC_many_swsperiods(Uscc, Dscc, 1)
    du = SCC_many_swsperiods(Uscc, Dscc, 0)

    return mean_up, cv_up, mean_down, cv_down, ud, du, perc_up, perc_down



# Main function that can be run from the terminal

def main():
    """Note: it may be better to just use the calculate stats function above
    
    Wrapper function that takes in the data and prints out the statistics (mean up, mean down, cv up, cv down, scc u->d, scc d->u) for each recording in the data
    
    sys.argv contains:
        up_file   -- file name of a csv file containing start and stop times for many up states
        down_file -- file name of a csv file containing start and stop times for many down states
        sws_file  -- file name of a csv file containing start and stop times of slow wave sleep episodes
            previously, I used these instead of sws_file:
            sws_start -- start time of sws period
            sws_stop  -- stop time of sws period
    
    Can also have sys.argv contain just up_file and down_file if uninterested in SCCs
    """
    if len(sys.argv) == 3:
        # Process input
        up_file = sys.argv[1]
        down_file = sys.argv[2]


        # Create lists of up and down durations

        up_times = list(csv.reader(open(up_file)))
        down_times = list(csv.reader(open(down_file)))

        up_dur = []
        down_dur = []
        for timePairU in up_times[1:]:
            state_start = float(timePairU[0])
            state_stop = float(timePairU[1])
            dur = state_stop - state_start
            up_dur.append(dur)
        for timePairD in down_times[1:]:
            state_start = float(timePairD[0])
            state_stop = float(timePairD[1])
            dur = state_stop - state_start
            down_dur.append(dur)

        # Calculate and print out statistics
        mean_up, mean_down, cv_up, cv_down, skew_up, skew_down = mean_cv_skew(up_dur, down_dur)
        scc = SCC(up_dur, down_dur, 1)
        print("Mean up: ", mean_up)
        print("CV up: ", mean_down)
        print("Mean down: ", cv_up)
        print("CV down: ", cv_down)
        #print("SCC: ", scc)
        #print("Skew up: ", skew_up)
        #print("Skew down: ", skew_down)

   
    elif len(sys.argv) == 4:
        # Process input
        up_file = sys.argv[1]
        down_file = sys.argv[2]
        sws_file = sys.argv[3]
        #sws_start = int(sys.argv[3])
        #sws_stop = int(sys.argv[4])

        # Create lists of up and down durations

        up_times_list = list(csv.reader(open(up_file)))
        down_times_list = list(csv.reader(open(down_file)))
        sws_times_list = list(csv.reader(open(sws_file)))

        up_times = np.array(up_times_list)
        down_times = np.array(down_times_list)
        sws_times = np.array(sws_times_list)
        
        # data for plotting:
        du_sccs = []
        ud_sccs = []
        num_ups = []

        Uscc = []
        Dscc = []
        Umcv = []
        Dmcv = []
        for swEpisode in sws_times[1:]:
            sws_start = float(swEpisode[0])
            sws_stop = float(swEpisode[1])

            up_dur = []
            down_dur = []
            for timePair in up_times[1:]:
                state_start = float(timePair[0])
                state_stop = float(timePair[1])
                if state_start >= sws_start and state_stop <= sws_stop:
                    dur = state_stop - state_start
                    up_dur.append(dur)
            for timePair in down_times[1:]:
                state_start = float(timePair[0])
                state_stop = float(timePair[1])
                if state_start >= sws_start and state_stop <= sws_stop:
                    dur = state_stop - state_start
                    down_dur.append(dur)
        
            Uscc.append(up_dur)
            Dscc.append(down_dur[:-1])

            Umcv.extend(up_dur)
            Dmcv.extend(down_dur)

        # Calculate and print out statistics
        mean_up, cv_up, mean_down, cv_down, skew_up, skew_down = mean_cv_skew(Umcv, Dmcv)
        ud = SCC_many_swsperiods(Uscc, Dscc, 1)
        du = SCC_many_swsperiods(Uscc, Dscc, 0)

        print("Mean up: ", mean_up)
        print("CV up: ", cv_up)
        print("Mean down: ", mean_down)
        print("CV down: ", cv_down)
        print("U->D SCC: ", ud)
        print("D->U SCC: ", du)



            # for plotting:
            #if math.isnan(ud) == False & math.isnan(du) == False:
            #    du_sccs.append(du)
            #    ud_sccs.append(ud)
            #    num_ups.append(len(up_dur_array))

        #print(du_sccs)
        #print(ud_sccs)
        #print(num_ups)
        #return(du_sccs, ud_sccs)

        #rows = zip(du_sccs, ud_sccs, num_ups)
        #with open("C:/Users/Shoshana/Documents/CSHL Summer/sws_stats_plots/scc_BWRat17_121912_real.csv", "w") as f:
        #    writer = csv.writer(f)
        #    for row in rows:
        #        writer.writerow(row)


    else:
        raise Exception("Incorrect number of arguments: %s" % sys.argv)

if __name__ == '__main__':
    sys.exit(main())