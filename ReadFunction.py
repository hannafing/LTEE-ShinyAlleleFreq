import pandas as pd
import numpy as np


# m1 --------------------------------------------------------------------------------------------
#Read the data from the file
initial_data_m1 = pd.read_csv("m1_annotated_timecourse.txt", sep=',')

initial_data_m1.to_csv('initial_data.csv', index=False)

#Mask data to only include rows where the 'Passed?' column is not 'FAIL'
initial_data_m1.columns = initial_data_m1.columns.str.strip()
initial_data_m1['Passed?'] = initial_data_m1['Passed?'].str.strip()
initial_data_m1 = initial_data_m1[initial_data_m1['Passed?'] != 'FAIL']

#print(initial_data.head())

data_clean_m1 = initial_data_m1.copy()

#Save the cleaned data to a new CSV file
data_clean_m1.to_csv("data_clean.csv", index=False)


#Create final dataframe
data_m1 = data_clean_m1[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

#Turn all column names into list / initiate empty list for AC columns
all_columns = list(data_clean_m1.columns)
ac_cols = []

#Loop through all columns and append ones that start with 'AC:'
for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

#Loop through each AC column and calculate allele frequency
for ac_col in ac_cols:

    #Split column name into two parts
    parts = ac_col.split(':')
    timepoint = parts[1]

    #Make matching column names for DP and AF
    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    #Replace any 0 depth values with NaN (no division by 0))
    depth = data_clean_m1[dp_col].replace(0, float('nan'))

    #Calculate allele frequency and store it as a new column
    data_m1[af_col] = data_clean_m1[ac_col] / depth

data_m1.to_csv("data_m1.csv", index=False)


# m2 --------------------------------------------------------------------------------------------
initial_data_m2 = pd.read_csv("m2_annotated_timecourse.txt", sep=',')

initial_data_m2.to_csv('initial_data.csv', index=False)

initial_data_m2.columns = initial_data_m2.columns.str.strip()
initial_data_m2['Passed?'] = initial_data_m2['Passed?'].str.strip()
initial_data_m2 = initial_data_m2[initial_data_m2['Passed?'] != 'FAIL']

data_clean_m2 = initial_data_m2.copy()

data_clean_m2.to_csv("data_clean.csv", index=False)


data_m2 = data_clean_m2[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_m2.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)


for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_m2[dp_col].replace(0, float('nan'))

    data_m2[af_col] = data_clean_m2[ac_col] / depth

data_m2.to_csv("data_m2.csv", index=False)


# m3 --------------------------------------------------------------------------------------------
initial_data_m3 = pd.read_csv("m3_annotated_timecourse.txt", sep=',')

initial_data_m3.to_csv('initial_data.csv', index=False)

initial_data_m3.columns = initial_data_m3.columns.str.strip()
initial_data_m3['Passed?'] = initial_data_m3['Passed?'].str.strip()
initial_data_m3 = initial_data_m3[initial_data_m3['Passed?'] != 'FAIL']

data_clean_m3 = initial_data_m3.copy()

data_clean_m3.to_csv("data_clean.csv", index=False)


data_m3 = data_clean_m3[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_m3.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_m3[dp_col].replace(0, float('nan'))

    data_m3[af_col] = data_clean_m3[ac_col] / depth

data_m3.to_csv("data_m3.csv", index=False)


# m4 --------------------------------------------------------------------------------------------
initial_data_m4 = pd.read_csv("m4_annotated_timecourse.txt", sep=',')

initial_data_m4.to_csv('initial_data.csv', index=False)

initial_data_m4.columns = initial_data_m4.columns.str.strip()
initial_data_m4['Passed?'] = initial_data_m4['Passed?'].str.strip()
initial_data_m4 = initial_data_m4[initial_data_m4['Passed?'] != 'FAIL']

data_clean_m4 = initial_data_m4.copy()

data_clean_m4.to_csv("data_clean.csv", index=False)

data_m4 = data_clean_m4[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_m4.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_m4[dp_col].replace(0, float('nan'))

    data_m4[af_col] = data_clean_m4[ac_col] / depth

data_m4.to_csv("data_m4.csv", index=False)


# m5 --------------------------------------------------------------------------------------------
initial_data_m5 = pd.read_csv("m5_annotated_timecourse.txt", sep=',')

initial_data_m5.to_csv('initial_data.csv', index=False)

initial_data_m5.columns = initial_data_m5.columns.str.strip()
initial_data_m5['Passed?'] = initial_data_m5['Passed?'].str.strip()
initial_data_m5 = initial_data_m5[initial_data_m5['Passed?'] != 'FAIL']

data_clean_m5 = initial_data_m5.copy()

data_clean_m5.to_csv("data_clean.csv", index=False)

data_m5 = data_clean_m5[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_m5.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_m5[dp_col].replace(0, float('nan'))

    data_m5[af_col] = data_clean_m5[ac_col] / depth

data_m5.to_csv("data_m5.csv", index=False)


# m6 --------------------------------------------------------------------------------------------
initial_data_m6 = pd.read_csv("m6_annotated_timecourse.txt", sep=',')

initial_data_m6.to_csv('initial_data.csv', index=False)

initial_data_m6.columns = initial_data_m6.columns.str.strip()
initial_data_m6['Passed?'] = initial_data_m6['Passed?'].str.strip()
initial_data_m6 = initial_data_m6[initial_data_m6['Passed?'] != 'FAIL']

data_clean_m6 = initial_data_m6.copy()

data_clean_m6.to_csv("data_clean.csv", index=False)

data_m6 = data_clean_m6[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_m6.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_m6[dp_col].replace(0, float('nan'))

    data_m6[af_col] = data_clean_m6[ac_col] / depth

data_m6.to_csv("data_m6.csv", index=False)


# p1 --------------------------------------------------------------------------------------------
initial_data_p1 = pd.read_csv("p1_annotated_timecourse.txt", sep=',')

initial_data_p1.to_csv('initial_data.csv', index=False)

initial_data_p1.columns = initial_data_p1.columns.str.strip()
initial_data_p1['Passed?'] = initial_data_p1['Passed?'].str.strip()
initial_data_p1 = initial_data_p1[initial_data_p1['Passed?'] != 'FAIL']

data_clean_p1 = initial_data_p1.copy()

data_clean_p1.to_csv("data_clean.csv", index=False)

data_p1 = data_clean_p1[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_p1.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_p1[dp_col].replace(0, float('nan'))

    data_p1[af_col] = data_clean_p1[ac_col] / depth

data_p1.to_csv("data_p1.csv", index=False)


# p2 --------------------------------------------------------------------------------------------
initial_data_p2 = pd.read_csv("p2_annotated_timecourse.txt", sep=',')

initial_data_p2.to_csv('initial_data.csv', index=False)

initial_data_p2.columns = initial_data_p2.columns.str.strip()
initial_data_p2['Passed?'] = initial_data_p2['Passed?'].str.strip()
initial_data_p2 = initial_data_p2[initial_data_p2['Passed?'] != 'FAIL']

data_clean_p2 = initial_data_p2.copy()

data_clean_p2.to_csv("data_clean.csv", index=False)

data_p2 = data_clean_p2[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_p2.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_p2[dp_col].replace(0, float('nan'))

    data_p2[af_col] = data_clean_p2[ac_col] / depth

data_p2.to_csv("data_p2.csv", index=False)


# p3 --------------------------------------------------------------------------------------------
initial_data_p3 = pd.read_csv("p3_annotated_timecourse.txt", sep=',')

initial_data_p3.to_csv('initial_data.csv', index=False)

initial_data_p3.columns = initial_data_p3.columns.str.strip()
initial_data_p3['Passed?'] = initial_data_p3['Passed?'].str.strip()
initial_data_p3 = initial_data_p3[initial_data_p3['Passed?'] != 'FAIL']

data_clean_p3 = initial_data_p3.copy()

data_clean_p3.to_csv("data_clean.csv", index=False)

data_p3 = data_clean_p3[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_p3.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_p3[dp_col].replace(0, float('nan'))

    data_p3[af_col] = data_clean_p3[ac_col] / depth

data_p3.to_csv("data_p3.csv", index=False)


# p4 --------------------------------------------------------------------------------------------
initial_data_p4 = pd.read_csv("p4_annotated_timecourse.txt", sep=',')

initial_data_p4.to_csv('initial_data.csv', index=False)

initial_data_p4.columns = initial_data_p4.columns.str.strip()
initial_data_p4['Passed?'] = initial_data_p4['Passed?'].str.strip()
initial_data_p4 = initial_data_p4[initial_data_p4['Passed?'] != 'FAIL']

data_clean_p4 = initial_data_p4.copy()

data_clean_p4.to_csv("data_clean.csv", index=False)

data_p4 = data_clean_p4[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_p4.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_p4[dp_col].replace(0, float('nan'))

    data_p4[af_col] = data_clean_p4[ac_col] / depth

data_p4.to_csv("data_p4.csv", index=False)


# p5 --------------------------------------------------------------------------------------------
initial_data_p5 = pd.read_csv("p5_annotated_timecourse.txt", sep=',')

initial_data_p5.to_csv('initial_data.csv', index=False)

initial_data_p5.columns = initial_data_p5.columns.str.strip()
initial_data_p5['Passed?'] = initial_data_p5['Passed?'].str.strip()
initial_data_p5 = initial_data_p5[initial_data_p5['Passed?'] != 'FAIL']

data_clean_p5 = initial_data_p5.copy()

data_clean_p5.to_csv("data_clean.csv", index=False)

data_p5 = data_clean_p5[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_p5.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_p5[dp_col].replace(0, float('nan'))

    data_p5[af_col] = data_clean_p5[ac_col] / depth

data_p5.to_csv("data_p5.csv", index=False)


# p6 --------------------------------------------------------------------------------------------
initial_data_p6 = pd.read_csv("p6_annotated_timecourse.txt", sep=',')

initial_data_p6.to_csv('initial_data.csv', index=False)

initial_data_p6.columns = initial_data_p6.columns.str.strip()
initial_data_p6['Passed?'] = initial_data_p6['Passed?'].str.strip()
initial_data_p6 = initial_data_p6[initial_data_p6['Passed?'] != 'FAIL']

data_clean_p6 = initial_data_p6.copy()

data_clean_p6.to_csv("data_clean.csv", index=False)

data_p6 = data_clean_p6[['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 
                   'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 
                   'Fold increase', 'Duplication pvalue', 'Passed?']].copy()

all_columns = list(data_clean_p6.columns)
ac_cols = []

for col in all_columns:
    if col.startswith('AC:'):
        ac_cols.append(col)

for ac_col in ac_cols:

    parts = ac_col.split(':')
    timepoint = parts[1]

    dp_col = 'DP:' + timepoint
    af_col = 'AF:' + timepoint

    depth = data_clean_p6[dp_col].replace(0, float('nan'))

    data_p6[af_col] = data_clean_p6[ac_col] / depth

data_p6.to_csv("data_p6.csv", index=False)