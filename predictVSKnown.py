from datetime import datetime, timedelta

# The epidemiological information for each patient.
class Patient(object):

    def __init__(self, p_id, p_inf_d, p_inf_id): # Patient initializer
        self.id = p_id
        self.cluster = ""
        self.p_inf_d = p_inf_d
        self.lntd = None
        self.fptd = None
        self.p_inf_id = p_inf_id
        self.sex = ""
        self.clade = ""
        self.a_inf_id = ""
        self.TEv = 0
        self.TType = 0
        self.IEv = 0
        self.IType = 0

    def add_BaseInfo(self, cluster, lntd, fptd, sex, clade, a_inf_id): 
        self.cluster = cluster
        self.lntd = lntd
        self.fptd = fptd
        self.sex = sex
        self.clade = clade
        self.a_inf_id = a_inf_id

    def report(self, pt_array, pt_dict):
        if (self.a_inf_id == "NA" or self.a_inf_id == ""):
            info = str(self.id) + "\t" + str(self.cluster) + "\t" + str(self.clade) + "\t" + str(self.sex) + "\t" + str(self.TEv) + "\t" + str(self.TType) + "\t" + str(self.IEv) + "\t" + str(self.IType) + "\t" + str(self.a_inf_id) + "\tNA\t" + str(self.p_inf_id) + "\t" + str(self.lntd) + "\t" + str(self.fptd) + "\t" + str(self.p_inf_d) + "\n"
            print(str(self.id) + "\t" + str(self.p_inf_d))
        else:
            tmp_ind = pt_dict[self.a_inf_id]
            # print(self.a_inf_id + "\t" + str(tmp_ind))
            info = str(self.id) + "\t" + str(self.cluster) + "\t" + str(self.clade) + "\t" + str(self.sex) + "\t" + str(self.TEv) + "\t" + str(self.TType) + "\t" + str(self.IEv) + "\t" + str(self.IType) + "\t" + str(self.a_inf_id) + "\t" + str(pt_array[int(tmp_ind)].sex) + "\t" + str(self.p_inf_id) + "\t" + str(self.lntd) + "\t" + str(self.fptd) + "\t" + str(self.p_inf_d) + "\n"
            print(str(self.id) + "\t" + str(self.p_inf_d))
        return info

def add_TestDate(entry):
    if (entry == "NA"): # Have preset all empty cells to NA
        return None
    else:
        return datetime.strptime(entry, '%m/%d/%Y').date()

wkdir = "./" # "/home/lmak/Dropbox/University of Calgary/SAC HIV Deep Sequencing/Transmission Sequences/"

# 1) Load the inferred transmissions and store all of the information. 
file = wkdir + "predicted_transmission_table.tsv"
inferred = open(file)
line = inferred.readline() # Skip the header row. 
line = inferred.readline()
pt_array = []
pt_dict = {} # pt_dict and ct match the patient ID (key) to the array index (value). 
ct = 0
while (ct < 140):
    p_info = line.split('\t')
    p_id = p_info[1].split('_')[0] # 5459_10-Jan-18_Prrt_14262
    p_pinfd = datetime.strptime(p_info[2].strip(), '%Y-%m-%d').date() # 2016-07-11
    p_pinfid = "NA"
    if (len(p_info[4]) > 3):
        p_pinfid = p_info[4].split('_')[0].strip() # 5462_12-Jan-18_Prrt_14296
    pt_array.append(Patient(p_id, p_pinfd, p_pinfid))
    pt_dict[p_id] = ct
    ct += 1 
    line = inferred.readline()
inferred.close()

# 2) Load the actual transmissions and store all of the relevant information. 
file = wkdir + "actual_transmission_table.tsv"
actual = open(file)
line = actual.readline() # Skip the header row. 
line = actual.readline()
cl_dict = {} # cl_dict keeps track of the clusters and will be used to evaluate actual transmission events. 
prev_id = ""
while line: # aa    Patient_ID    14-05-2002    06688A    I    01-01-1992    01-01-1993    B
    p_info = line.split('\t')
    p_id = p_info[1]
    if (p_id is prev_id):
        continue
    prev_id = p_id    
    p_cluster = p_info[0].strip()
    if (p_cluster not in cl_dict): 
        cl_dict[p_cluster] = [] 
    cl_dict[p_cluster].append(p_id) # Makes looking for actual transmissions easier.
    p_sex = p_info[4]
    p_lntd = add_TestDate(p_info[5].strip())
    p_fptd = add_TestDate(p_info[6].strip())
    p_clade = p_info[7]
    p_ainfid = p_info[8]
    pt_array[pt_dict[p_id]].add_BaseInfo(p_cluster, p_lntd, p_fptd, p_sex, p_clade, p_ainfid)
    line = actual.readline()
actual.close() # All of the patient objects are now finished.

# 3) Accuracy of transmission and infection date inferences? Also, write to a report file.
file = wkdir + "transmission_report.txt"
report = open(file,'w')
report.write("SAC ID\tCluster\tClade\tSex\tTrans. Ev.\tTrans. Acc.\tInfect. Ev.\tInfect. Acc.\tActual Infector\tSex\tPred. Infector\tLast Negative Date\tFirst Positive Date\tPred. Infect. Date\n")
for pt in pt_array:
    if (pt.cluster != "NA"): 
        p_inf_id_cl = ""
        try:
            p_inf_id_cl = pt_array[pt_dict[pt.p_inf_id]].cluster
        except:
            print()
        if (pt.a_inf_id != "NA"): # 1) Transmission evidence? Has i) cluster and ii) is not the originator.
            pt.TEv = 2
            if (pt.p_inf_id == pt.a_inf_id): # 1. Right: Correct ID. Will have to manually check for uncertain cases. Patient.P_Inf_ID == Patient.A_Inf_ID. 
                pt.TType = 1
            else:
                if (pt.cluster == p_inf_id_cl): # 2. Close: Correct cluster, wrong direction. Above not true but Patient.Cluster == Patient.P_Inf_ID.Cluster.
                    pt.TType = 2
                else: # 0. Wrong: Neither of the above. 
                    pt.TType = 0
        else: # 2) No transmission evidence because is cluster originator.  
            # print(str(pt.id) + "\t" + pt.cluster + "\t" + p_inf_id_cl)
            pt.TEv = 1
            if (pt.cluster == p_inf_id_cl): # 2. Close: They are i) the originator of a cluster and ii) infected by someone in the same cluster. Patient.Cluster == Patient.P_Inf_ID.Cluster.
                pt.TType = 2
            else: # 3. Totally new inference. May be a patient from another cluster or a non-clustered patient.
                pt.TType = 3
    else: # 2) No transmission evidence because is non-clustered. 
        pt.TEv = 0
        pt.TType = 3
    # print(str(pt.id) + "\t" + str(pt.lntd) + "\t" + str(pt.fptd))
    if (pt.lntd is not None):
        pt.IEv = 2
        if (pt.fptd is not None): # 3) If lntd and fptd datetime.date() objects are not empty...
            if (pt.lntd <= pt.p_inf_d <= pt.fptd): # 1. Right: LNTD <= P_Inf_D <= FPTD.
                pt.IType = 1 
            else: # 0. Wrong: Above not true.
                pt.IType = 0
    else: 
        if (pt.fptd is not None): # 4) If only lntd is empty:
            pt.IEv = 1
            if (pt.p_inf_d <= pt.fptd): # 3. New: P_Inf_D <= FPTD
                pt.IType = 1
            else: # Wrong: Above not true. 
                pt.IType = 0
        else: # 5) No information. 
            pt.IEv = 0
            pt.IType = 3        
    report.write(pt.report(pt_array, pt_dict))
report.close()