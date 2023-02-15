import matplotlib.pyplot as plt
import pandas as pd

num_list=[]
for i in range(2,51):
    print(i)
    data_df = pd.read_csv("../results/mufuzz_consecutive/cluster/"+ str(i) + "/go_number.txt",sep="\t")
    number = data_df.iloc[i,0]
    num_list.append(number)

plt.plot(num_list)
plt.xlabel("Number of clusers")
plt.ylabel("Total number of ORA GO teams")