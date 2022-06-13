import csv
import pandas as pd
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

# dataset = pd.read_csv('play_tennis.csv')
# output_file = "onehot_tennis.dat"
dataset = pd.read_csv('kohkiloyeh.csv')
output_file = "onehot_kohkiloyeh.dat"
x = dataset.iloc[:, :].values
print(x)

labelencoder = LabelEncoder()
lenlist = []
for i in range(len(x[0])):
    x[:, i] = labelencoder.fit_transform(x[:, i])
    lenlist.append(len(set(x[:, i])))
    print(lenlist[i])
onehot_encoder = OneHotEncoder(sparse=False, dtype=int)

result = []
bitlen = 0
with open(output_file, "w") as out:
    # writer = csv.writer(out)
    for i in range(len(x[0])):
        temp_cv = x[:, i].reshape(len(x[:, i]), 1)
        if lenlist[i] > 2:
            bitlen += lenlist[i]
            temp_binary = onehot_encoder.fit_transform(temp_cv)
            result.append(temp_binary)
        else:
            bitlen += 1
            result.append(temp_cv)
    out.write('%d %d %d\n' % (len(result), bitlen, len(result[0])))
    for i in range(len(result[0])):
        strA = ""
        for j in range(len(result)):
            str_temp = " ".join([str(x) for x in result[j][i]])
            # out.write('%d ' % str_temp)
            strA += str_temp
            strA += " "
        strA += "\n"
        out.writelines(strA)
        print(strA)
