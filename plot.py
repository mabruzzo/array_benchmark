import glob2
import csv
import numpy as np
import matplotlib.pyplot as plt

use_errs = False
use_min = False

results = []
N_vals = []
for fname in sorted(list(glob2.glob("result_*.txt"))):
    N_vals.append( int(fname[7:-4]))
    results.append(np.loadtxt(fname, delimiter=", "))
combined = np.stack(results)

test_names = ["pointer", "boost::multi_array[][][]",
              "boost::multi_array()",
              "EnzoArray", "EnzoArray_neg"]
data = {}
N_vals = np.array(N_vals)
idx = np.argsort(N_vals)
fmts = ["ko:","rs","yo","g>", "b."]

for i,name in enumerate(test_names):
    cur_data = combined[:,i,:]
    if use_errs:
        cur_vals = np.mean(cur_data,axis=1)
        cur_err = np.std(cur_data,axis=1)
    else:
        cur_vals = np.amin(cur_data,axis=1)
        cur_err = np.zeros_like(cur_vals)
    data[name] = np.hstack([cur_vals[:,np.newaxis],cur_err[:,np.newaxis]])


fig,ax = plt.subplots(2,1, figsize=(4,7),sharex=True)
for name,fmt in zip(test_names,fmts):
    ax[0].errorbar(N_vals[idx],data[name][:,0][idx],
                   yerr=data[name][:,1][idx], fmt=fmt,label = name)
ax[0].set_yscale("log")
ax[0].set_xscale("log",basex=2)
ax[0].set_ylabel("$t(N^3)$")
ax[0].legend()

for name,fmt in zip(test_names,fmts)[1:]:
    ratio = data[name][:,0]/data["pointer"][:,0]
    err = ratio*np.sqrt(np.square(data[name][:,1]/data[name][:,0]) +
                         np.square(data["pointer"][:,1]/data["pointer"][:,0]))
    ax[1].errorbar(N_vals[idx], ratio[idx], yerr=err[idx], fmt=fmt,
                   label = name, capsize=2.0)
ax[1].set_ylabel(r"$t(N^3)/t_{\rm pointer}(N^3)$")
ax[1].set_xlabel("$N$")
ax[1].axhline(1.0, color = "grey", linestyle = "--")
fig.tight_layout()
plt.savefig("plot.pdf")
