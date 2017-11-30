import random as rd
import csv

def main():
    # Attention gatherers
    AGs = ["AG1", "AG2", "AG3", "AG4", "AG5"]
    # Category names
    tail_cat = ["A", "B"]
    head_cat = ["1", "2"]
    # t/csv file headers
    header_cb = ["Participant"] + ["Fam"+str(i) for i in range(8)] + ["AG"]
    header_cb += ["Fam"+str(i) for i in range(8,16)] + ["AG"]
    header_cb += ["Fam"+str(i) for i in range(16,24)] + ["AG"]
    header_cb += ["Ctr0", "AG", "Ctr1", "AG", "Ctr2", "AG", "WL0", "AG", "WL1"]
    header_pi = ["Participant", "First fam stim", "Name first fam stim",
                 "First contrast", "Side new HC", "Side new TC", "Side new tail RC",
                 "First label WL", "Side A1", "Side A2"]
    rows_cb = [header_cb]
    rows_pi = [header_pi]
    # Label condition (twice NoLabel to get even numbers
    # in both Label and NoLabel conditions)
    conditions = {"Label":["Gatoo", "Saldie"],
                  "NoLabel":["NoLabel", "NoLabel"]
                 }
    # Familiarisation stimlus names
    fam_stims = [{"Gatoo":[["A"+str(i)+"_G1" for i in range(1,7)],
                           ["A"+str(i)+"_G2" for i in range(1,7)]],
                  "Saldie":[["A"+str(i)+"_S1" for i in range(1,7)],
                            ["A"+str(i)+"_S2" for i in range(1,7)]],
                  "NoLabel":[["A"+str(i)+"_NL1" for i in range(1,7)],
                             ["A"+str(i)+"_NL2" for i in range(1,7)]]
                 },
                 {"Gatoo":[["B"+str(i)+"_G1" for i in range(1,7)],
                           ["B"+str(i)+"_G2" for i in range(1,7)]],
                  "Saldie":[["B"+str(i)+"_S1" for i in range(1,7)],
                            ["B"+str(i)+"_S2" for i in range(1,7)]],
                  "NoLabel":[["B"+str(i)+"_NL1" for i in range(1,7)],
                             ["B"+str(i)+"_NL2" for i in range(1,7)]]
                 }
                ]
    # Familiarisation stimulus order
    order1 = list(range(6))
    order2 = list(range(6))
    # Familiarisation stimulus facing direction
    dirs = ["Reg_", "Flip_"]
    # Familiarisation counterbalancing
    p = 0
    while p < 48:
        for condition in conditions:
            for first_stim in (0,1):
                for name_first_stim in (0,1):
                    # Generte random order for each category
                    rd.shuffle(order1)
                    rd.shuffle(order2)
                    # Increment participant id (to start with participant 1)
                    p += 1
                    row_cb = [p]
                    row_pi = [p]
                    for i in range(12):
                        # Repeat the order twice to get 12 pres/cat in total, with AGs
                        if i %4 == 0 and i > 0:
                            row_cb.append(AGs[2 + i//4])
                        row_cb.append(rd.choice(dirs)+fam_stims[first_stim][conditions[condition][name_first_stim]][i//6][order1[i%6]])
                        row_cb.append(rd.choice(dirs)+fam_stims[first_stim-1][conditions[condition][name_first_stim-1]][i//6][order2[i%6]])
                    row_pi.append([tail_cat[first_stim], conditions[condition][name_first_stim]])
                    rows_cb.append(row_cb)
                    rows_pi.append(row_pi)
    # Contrast tests
    cb_tests = ["HC", "TC"]
    Oh_HC_sides = ["L","R"]
    Ot_TC_sides = ["L","R"]
    Ot_RC_sides = ["L","R"]
    # Contrast tests counterbalancing
    p = 0
    while p < 48:
        for first_test in (0,1):
            for Oh_HC in (0,1):
                for Ot_TC in (0,1):
                    for Ot_RC_side in (0,1):
                        p += 1
                        Ot_HC = rd.choice([0,1])
                        Oh_HC_side = rd.choice([0,1])
                        Oh_TC = rd.choice([0,1])
                        Ot_TC_side = rd.choice([0,1])
                        HC_TC = {"HC":"HC_"+tail_cat[Ot_HC]+head_cat[Oh_HC]+Oh_HC_sides[Oh_HC_side],
                                 "TC":"TC_"+tail_cat[Ot_TC]+head_cat[Oh_TC]+Ot_TC_sides[Ot_TC_side]
                                }
                        RC = "RC_"+tail_cat[Ot_TC-1]+Ot_RC_sides[Ot_RC_side]+"_"+head_cat[Oh_HC-1]+Ot_RC_sides[Ot_RC_side-1]
                        rows_cb[p] += [AGs[0], HC_TC[cb_tests[first_test]], AGs[1], HC_TC[cb_tests[first_test-1]], AGs[2], RC]
                        rows_pi[p] += [cb_tests[first_test], Oh_HC_sides[Oh_HC_side-1],
                                       Ot_TC_sides[Ot_TC_side-1], Ot_RC_sides[Ot_RC_side-1]]
    
    # Word-learning tests
    labels = ["G","S"]
    A1_sides = ["L","R"]
    A2_sides = ["L","R"]
    # Word learning tests counterbalancing
    p = 0
    while p < 48:
        for first_label in (0,1):
            for A1_side in (0,1):
                for A2_side in (0,1):
                    p += 1
                    rows_cb[p] += [AGs[3], "WL"+labels[first_label]+"_A1"+A1_sides[A1_side]+"_B2"+A1_sides[A1_side-1],
                                   AGs[4], "WL"+labels[first_label-1]+"_A2"+A2_sides[A2_side]+"_B1"+A2_sides[A2_side-1]]
                    rows_pi[p] += [labels[first_label], A1_sides[A1_side], A2_sides[A2_side]]
    
    with open('Counterbalancing.tsv', 'w') as cb, open('ParticipantInfo.csv', 'w') as pi:
        cb = csv.writer(cb, delimiter='\t')
        pi = csv.writer(pi)
        cb.writerows(rows_cb)
        pi.writerows(rows_pi)

if __name__=="__main__":
    main()
