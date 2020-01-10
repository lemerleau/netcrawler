import numpy as np 
import pandas 
import matplotlib.pyplot as plt 



def main(): 
	"""
	df1 = pandas.read_csv("popbase_data.csv")
	df2 = pandas.read_csv("netcrawlerbase_data.csv")
	plt.rcParams['axes.labelweight'] = 'bold'
	plt.title("Average max fitness over 100 runs", fontname="Times New Roman Bold")
	plt.plot(df1['1'], label="population base")
	plt.plot(df2['1'], label="netcrawler")
	plt.legend(loc="lower right")
	plt.xlabel("Generation")
	plt.ylabel("Max Fitness")
	plt.savefig("netcrawlerVSpop.eps")
	plt.savefig("netcrawlerVSpop.png")
	plt.show()
	"""
	df_all = pandas.read_csv("all_pop_data.txt")
	#df_all.rename(columns={'0':"id",'0': "Generation", '1': "popBase"})

	df_all['1'].plot.hist(bins=12, label='Pop base')
        
	df_all2 = pandas.read_csv("all_netcrawler_data.txt")
	df_all2['1'].plot.hist(bins=14, label='netcrawler')
	plt.legend(loc='upper right')
	plt.title("Max fitness Distribution")
	plt.savefig("full_data_hist.png")
	plt.xlabel("Max fitness")
	plt.show()




if __name__ =="__main__": 

	main()
