SEQUENCE_1 = 'GCCAAGACAGTATCAGACAC' #Sequence 1 (Side Sequence)
SEQUENCE_2 = 'TGTCATGCCAAGACAGTATG' #Sequence 2 (Top Sequence)
MATRIX_ROW_N = len(SEQUENCE_1)+1 #Initiation Matrix Size (Rows)
MATRIX_COLUMN_N = len(SEQUENCE_2)+1 #Initiation Matrix Size (Columns)
MATCH_SCORE = 5 #Match Score
MISMATCH_SCORE = -3 #Mismatch Score
GAP_SCORE = -5 #Gap Points
GAP_CHARACTER = '-'#Character to Represent Gaps in Final Alignemnts
ALN_PATHWAYS = [] #Initiating List of Discovered aln Pathways
MATRIX = [[[[None] for i in range(2)] for i in range(MATRIX_COLUMN_N)] for i in range(MATRIX_ROW_N)] # Initiating Score Matrix

def print_aln_details(aln): # Function to print steps of particular aln
    print('\n>>>>>>>>>>>>>>>> Aln #:'+str(aln[4])+' <<<<<<<<<<<<<<<<\n\n------------------Steps------------------\n')
    for elem in aln[3]:
        print('Step: '+str(elem[0])+' -> Score = ' + str(elem[1])    + ' -> Dir: ' + ('\u2190' if elem[2]==1 else ('\u2196' if elem[2]==2 else ('\u2191' if elem[2]==3 else 'Error'))))
    print('\n-------------Final aln-------------\n')
    print(aln[1]+'\n'+aln[0])
    return

def print_all(ALIGNMENTS): #Function to print everything
    print('Total Alignments: ' + str(len(ALIGNMENTS)))
    print('Overall Score: '+str(ALIGNMENTS[0][3][0][1])+'\n')
    for elem in ALIGNMENTS:
        print_aln_details(elem)
    return
def print_alns_only(ALIGNMENTS): #Function to print only ALIGNMENTS
    print('Total Alignments: ' + str(len(ALIGNMENTS)))
    print('Overall Score: '+str(ALIGNMENTS[0][3][0][1])+'\n')
    for elem in ALIGNMENTS:
        print(elem[0]+'\n'+elem[1]+'\n')
    return
def find_each_path(c_i,c_j,path=''): #Nested function to discover new aln pathways
    global ALN_PATHWAYS 
    i = c_i 
    j = c_j 
    if i == 0 and j==0: 
        ALN_PATHWAYS.append(path) 
        return 2 
    dir_t = len(MATRIX[i][j][1]) 
    while dir_t<=1: 
        n_dir = MATRIX[i][j][1][0] if (i != 0 and j != 0) else (1 if i == 0 else (3 if j==0 else 0)) 
        path = path + str(n_dir) 
        if n_dir == 1: 
            j=j-1
        elif n_dir == 2:
            i=i-1
            j=j-1
        elif n_dir == 3:
            i=i-1
        dir_t = len(MATRIX[i][j][1])
        if i == 0 and j==0:
            ALN_PATHWAYS.append(path)
            return 3
    if dir_t>1:
        for dir_c in range(dir_t):
            n_dir = MATRIX[i][j][1][dir_c] if (i != 0 and j != 0) else (1 if i == 0 else (3 if j==0 else 0))
            tmp_path = path + str(n_dir)
            if n_dir == 1:
                n_i = i
                n_j=j-1
            elif n_dir == 2:
                n_i=i-1
                n_j=j-1
            elif n_dir == 3:
                n_i=i-1
                n_j = j
            find_each_path(n_i,n_j,tmp_path)
    return len(ALN_PATHWAYS)

#Main Code
#Matrix Evaulation [start]
for i in range(MATRIX_ROW_N):
    MATRIX[i][0] = [GAP_SCORE*i,[]]
for j in range(MATRIX_COLUMN_N):
     MATRIX[0][j] = [GAP_SCORE*j,[]]
for i in range(1,MATRIX_ROW_N):
    for j in range(1,MATRIX_COLUMN_N):
        score = MATCH_SCORE if (SEQUENCE_1[i-1] == SEQUENCE_2[j-1]) else MISMATCH_SCORE
        h_val = MATRIX[i][j-1][0] + GAP_SCORE
        d_val = MATRIX[i-1][j-1][0] + score
        v_val = MATRIX[i-1][j][0] + GAP_SCORE
        o_val = [h_val, d_val, v_val]
        MATRIX[i][j] = [max(o_val), [i+1 for i,v in enumerate(o_val) if v==max(o_val)]] # h = 1, d = 2, v = 3
#Matrix Evaulation [end]
OVERALL_SCORE = MATRIX[i][j][0]
score = OVERALL_SCORE
l_i = i
l_j = j
ALIGNMENTS = []
tot_aln = find_each_path(i,j)
aln_count = 0
#Compiling alignments based on discovered matrix pathways
for elem in ALN_PATHWAYS:
    i = l_i-1
    j = l_j-1
    side_aln = ''
    top_aln = ''
    step = 0
    aln_info = []
    for n_dir_c in range(len(elem)):
        n_dir = elem[n_dir_c]
        score = MATRIX[i+1][j+1][0]
        step = step + 1
        aln_info.append([step,score,n_dir])
        if n_dir == '2':
            side_aln = side_aln + SEQUENCE_1[i]
            top_aln = top_aln + SEQUENCE_2[j]
            i=i-1
            j=j-1
        elif n_dir == '1':
            side_aln = side_aln + GAP_CHARACTER
            top_aln = top_aln + SEQUENCE_2[j]
            j=j-1
        elif n_dir == '3':
            side_aln = side_aln + SEQUENCE_1[i]
            top_aln = top_aln + GAP_CHARACTER
            i=i-1
    aln_count = aln_count + 1
    ALIGNMENTS.append([top_aln[::-1],side_aln[::-1],elem,aln_info,aln_count])

print_alns_only(ALIGNMENTS)