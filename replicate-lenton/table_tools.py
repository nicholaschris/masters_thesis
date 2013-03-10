#!usr/bin/python
#table_tools.py

'''
Copy and paste from /home/nicholas/projects/skeleton/replicate-lenton/table_tools.py
NOTES:
Creates a table like this:

-----------------------------------------------------------------------------------------------------------------------------------------
|| Year      | Total simulated uptake | Sampling Estimate of Uptake | Sampling uncertainty ||
-----------------------------------------------------------------------------------------------------------------------------------------
|| year[0] | tot_sam_upt[0]            | sam_est_upt[0]                     | sam_unc[0]                ||
-----------------------------------------------------------------------------------------------------------------------------------------
||             |                                    |                                              |                                   ||
-----------------------------------------------------------------------------------------------------------------------------------------

'''

import os
import cPickle
import numpy as np
import sampling_tools
from numpy import ma

'''
import replicate_lenton_class
def make_a_new_lenton():
	new_lenton = replicate_lenton_class.Lenton_2006()
	new_lenton.make_array_from_nc_variable()
	new_lenton.year_stack_func()
	new_lenton.signal_and_noise_constructor()
	new_lenton.summer_and_winter_climatological()
	new_lenton.snr_assess_temporal_variability()
	new_lenton.make_2D_time_space_array()
	return new_lenton

#new_lenton  = make_a_new_lenton()
'''

filename = 'regridded_cflux.pkl'
#filename = 'cflux_5day.pkl'
cur_dir = os.getcwd()
f = open(cur_dir + '/../data/pkl_files/' + filename)
data = cPickle.load(f)
f.close()

year_stack = ma.array(np.split(data[:, :45, :], 10, axis=0))
# unit_conversation = 12 * area * 31536000 / 1e15 # mol * area * seconds / petagrams

### Get the years from the lenton_class
year_start = 1998
#year_start = 2000
years = range(year_start, year_start+10, 1)
time_end = 360
#time_end = 72
time_frequency = 72 
#time_frequency = 18
lon_frequency = 10 # 1 = 2degrees
sam_est_upt = []
sam_unc = []
sam_unc_2 = []

### sample the array in the lenton class to get the following lists.
#~ tot_sam_upt = range(10)
tot_sam_upt = np.round([sampling_tools.sample_all_realizations(year_stack[x, :time_end, :, :180], time_frequency=1, lon_frequency=1)[0] for x in range(10)], 2)
#~ sam_est_upt =  range(10)
sample = np.round([sampling_tools.sample_all_realizations(year_stack[x, :time_end, :, :180], time_frequency=time_frequency, lon_frequency=lon_frequency) for x in range(10)], 2)
for item in sample:
    sam_est_upt = np.append(sam_est_upt, item[0])
    sam_unc = np.append(sam_unc, 2*item[1]) # NB Multiplying by 2
    sam_unc_2 = np.append(sam_unc_2, item[1]) 

#sam_unc = np.round([sampling_tools.sample_all_realizations(year_stack[x, :360, :, :180], time_frequency=18, lon_frequency=15)[1] for x in range(10)], 2)
#filename = raw_input("Enter filename to write table to: ")

#1st column
mean_tot_sam_upt = str(np.mean(tot_sam_upt))
mean_tot_sam_upt_unc = 2*np.round(np.std(tot_sam_upt)/np.sqrt(10), 2)

#2nd column
mean_sam_est_upt = str(np.mean(sam_est_upt))
mean_sam_est_upt_unc = str(2*np.round(np.std(sam_est_upt)/np.sqrt(10), 2)) # NB Why multply by 2 AGAIN?

#3rd column
#sam_unc for each year
mean_sam_unc = np.round(np.mean(sam_unc), 2) #average of the 2 sigma
sampling_error = np.round(np.mean(sam_unc)/np.sqrt(10), 2)
sampling_error_2 = np.round(np.mean(sam_unc_2)*(2/np.sqrt(10)), 2)

tot_sam_error = np.sqrt(mean_tot_sam_upt_unc**2 + sampling_error**2)

mean_tot_sam_upt_unc = str(mean_sam_est_upt_unc)
sampling_error = str(np.round(sampling_error, 2))
sampling_error_2 = str(np.round(sampling_error_2, 2))

fobj = open(os.getcwd() + '/output/'+ filename + '-'+str(time_frequency)+'-'+str(lon_frequency)+'.tex', 'w')

prop_or_curr = "Proposed" # Change this depending on which sampling strategy is used...

text = \
" \
\\begin{center} \n \
\\begin{table}[h] \n \
\centering \n \
\caption{Comparison of the total Simulated Uptake With the Uptake from Our \
"+ prop_or_curr + "(Sampling every "+str(time_frequency)+" days and every "+str(2*lon_frequency)+" degrees in longitude) \
Sampling and the Sampling Error Introduced. The sampling error is "+sampling_error+" and the total sampling error is "+str(np.round(tot_sam_error, 2))+"} \n \
\\begin{tabular}[tbp]{@{}lp{3cm}p{3cm}p{3cm}@{}} \\toprule \n \
\scriptsize{Year} & \scriptsize{Total Simulated Uptake, PgC/yr} \n \
& \scriptsize{Sample estimate of Uptake, PgC/yr} \n \
& \scriptsize{Sampling uncertainty (2$\sigma$), PgC/yr} \\\ \n \
\hline \n \
"
#OR#
text_start = " \
\\begin{table}[h] \n \
\caption{Comparison of the total Simulated Uptake With the Uptake from Our \n \
"+ prop_or_curr + "(Sampling every "+str(time_frequency)+" days and every "+str(2*lon_frequency)+" degrees in longitude) \
Sampling and the Sampling Error Introduced. The total sampling error is "+str(tot_sam_error)+"} \n \
\\begin{tabular}[h]{|p{2.5cm}|p{2.5cm}|p{2.5cm}|p{2.5cm}|} \n \
\hline \n \
\scriptsize{Year} & \scriptsize{Total Simulated Uptake, PgC/yr} \n \
& \scriptsize{Sample estimate of Uptake, PgC/yr} \n \
& \scriptsize{Sampling uncertainty (2$\sigma$), PgC/yr} \\ \n \
\hline \n \
"

fobj.write(text)

for year in range(10):
	fobj.write("\scriptsize{" +str(years[year])+ "}    &   \scriptsize{"+str(tot_sam_upt[year])+"} &  \scriptsize{"+str(sam_est_upt[year])+"} &  \scriptsize{"+str(sam_unc[year])+"} \\\\")
	fobj.write("\n")

text_end ="\hline \n \
\scriptsize{" +str(years[0])+ "-" +str(years[9])+ "}    \
&   \scriptsize{"+mean_tot_sam_upt+'$\pm$'+mean_tot_sam_upt_unc+"} \
&  \scriptsize{"+mean_sam_est_upt+'$\pm$'+mean_sam_est_upt_unc+"} \
&   \scriptsize{"+str(mean_sam_unc)+"} \\\ \n \
\\bottomrule \n \
\end{tabular} \n \
\end{table} \n \
\end{center}"

fobj.write(text_end)
fobj.close()

