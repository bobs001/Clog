# covid-19 statistics
# https://www.worldometers.info/coronavirus

# North America
#########################

tag = "US"
# date = "03/21/2020"
# pos = 19774
# dead = 275
# rec = 0 # unknown
date = "03/23/2020"
pos = 35079
dead = 458
rec = 178

# tag = "Canada"
# date = "03/23/2020"
# pos = 1472
# dead = 20
# rec = 14

# tag = "California"
# # date = "03/20/2020"
# # pos = 1057
# # dead = 20
# # rec = 0 # unknown
# date = "03/23/2020"
# pos = 1802
# dead = 35
# rec = 6

# tag = "New York"
# # date = "03/20/2020"
# # pos = 7102
# # dead = 35
# # rec = 0 # unknown
# date = "03/23/2020"
# pos = 16900
# dead = 150
# rec = 0

# tag = "Washington"
# # date = "03/20/2020"
# # pos = 1376
# # dead = 74
# # rec = 0 # unknown
# date = "03/23/2020"
# pos = 1996
# dead = 95
# rec = 124

# tag = "New Mexico"
# date = "03/20/2020"
# pos = 363
# dead = 4
# rec = 0 # unknown

# Europe
#########################

# tag = "Italy"
# # date = "02/27/2020"
# # pos = 5883
# # dead = 233
# # rec = 0 # unknown
# # date = "03/08/2020"
# # pos = 7375
# # dead = 366
# # rec = 0 # unknown
# date = "03/23/2020"
# pos = 59138
# dead = 5476
# rec = 7024

# tag = "France"
# date = "03/23/2020"
# pos = 16689
# dead = 674
# rec = 2200

# tag = "Spain"
# # date = "03/22/2020"
# # pos = 28603
# # dead = 1756
# # rec = 0 # unknown
# date = "03/23/2020"
# pos = 33089
# dead = 2206 
# rec = 3355

# tag = "Germany"
# # date = "03/22/2020"
# # pos = 24873
# # dead = 94
# # rec =  0 # unknown
# date = "03/23/2020"
# pos = 27181
# dead = 113
# rec =  422

# tag = "UK"
# date = "03/23/2020"
# pos = 5683
# dead = 289
# rec = 135

# tag = "Switzerland"
# date = "03/23/2020"
# pos = 8547
# dead = 118
# rec = 131

#Other
#########################

# tag = "S. Korea"
# date = "03/23/2020"
# pos = 8961
# dead = 111
# rec = 3166

# tag = "China"
# date = "03/22/2020"
# pos = 81093
# dead = 3270
# rec = 72703

# tag = "Japan"
# date = "03/23/2020"
# pos = 1101
# dead = 41
# rec = 235

# World
#########################

# tag = "World"
# date = "03/23/2020"
# pos = 351083
# dead = 15337
# rec = 100608

# calculate
##########################
USPop = 331002650. # wickipedia
drate = float(dead)/float(pos)
CFR = float(dead)/float(dead + rec)

print "reg  :", tag
print "date :", date
print "pos  :", pos
print "dead :", dead
print "drate: {drate:.2f}%".format(drate=drate*100)
print "CFR  : {cfr:.2f}% (case fatality rate)".format(cfr=CFR*100)
print "US projection: {:.2f} Million".format(USPop*drate/1e6)

