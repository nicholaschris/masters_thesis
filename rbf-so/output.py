# -*- coding: utf-8 -*-

print "maximum error",max_error
print "Total sum square error",tot_error



for k in range(72/2):
    print "Time slice: "+str(2*k+1) + ": " + str(sum((sample_locations[:,0] == 2*k)==True)) + "\t   Time slice: "+str(2*k+2) + ": " + str(sum((sample_locations[:,0] == 2*k+1)==True))    
print "Time slice: "+str(73) + ": " + str(sum((sample_locations[:,0] == 72)==True))
#broken
#all_pos.append(pos[0])




data_interp = zeros((data.shape[0],data.shape[1],data.shape[2]))
if mod(i+5,50)==0:
    cmd = 'mkdir ./Data'
    os.system(cmd)
    cmd = 'mkdir ./Data/'+str(i+5) + "Points"
    os.system(cmd)
    
    if i+5 == 50:
        f = open('./Data/errors' ,'w')
    else:
        f = open('./Data/errors' ,'a')
    f.write("============================== \niteration"  +str(i+5) +" \nmaximum error " +str(max_error) + "\nTotal sum square error " +str(tot_error) + "\n")
    f.close()
    
    #NOTE: uncomment the next four lines if output should be ASCII characters
    f = open('./Data/' + str(i+5) + "Points/locationsASCII",'w')
    for lm in range(sample_locations.shape[0]):
        f.write(str(int(sample_locations[lm,0]))+ "\t" + str(int(sample_locations[lm,1])) + "\t" + str(int(sample_locations[lm,2]))  + "\n")
    f.close()
    
    
    f = open('./Data/' + str(i+5) + "Points/locations.pkl",'w')
    cPickle.dump(sample_locations,f)
    f.close()
    
    
    f = open('./Data/' + str(i+5) + "Points/interpolatedData.pkl",'w')
    data_interp = zeros((data.shape[0],data.shape[1],data.shape[2]))
    for i in range(s_out.shape[0]):
        t = locations[valid_locations[i],0]
        lon = locations[valid_locations[i],1]
        lat = locations[valid_locations[i],2]
        data_interp[t,lon,lat] = s_out[i]    
        
    cPickle.dump(data_interp,f)
    f.close()
    