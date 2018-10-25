import math
import matplotlib.pyplot as plt
from numpy.random import normal
from scipy.integrate import odeint

#Screening Vars
stRNA_bindings = normal(140., scale=10., size=10)
rf1_bindings = normal(34.4, scale=5., size=10)

stRNA_unbindings = normal(60.23, scale=5., size=10)
rf1_unbindings = normal(0.19, scale=0.05, size=10)

initial_codons = [200., 400., 1000., 2000., 3200.]

#Vars
initial_codon_stRNA = 0.
initial_codon_RF1 = 0.

def sys(y, t):
     codon = y[0]
     codon_stRNA = y[1]
     codon_RF1 = y[2]
     #the model equations
     d_codon_dt = (stRNA_unbind * codon_stRNA + rf1_unbind * codon_RF1) - ((stRNA_bind + rf1_bind) * codon)
     d_codon_stRNA_dt = (stRNA_bind * codon) - (stRNA_unbind * codon_stRNA)
     d_codon_RF1_dt = (rf1_bind * codon) - (rf1_unbind * codon_RF1)
     return [d_codon_dt, d_codon_stRNA_dt, d_codon_RF1_dt]

def solve(init_cond, time):
  result = odeint(sys, init_cond, time)
  codon_res = result[:, 0]
  codon_stRNA_res = result[:,1]
  codon_RF1_res = result[:, 2]

  ratio = [(codon_res[i], codon_stRNA_res[i], codon_RF1_res[i]) for i in range(len(codon_res))]

  prob = math.pow((ratio[len(ratio) - 1][1]/initial_codon), (float(initial_codon/200.)))
  return prob, ratio

def run(initial_codon, initial_codon_stRNA, initial_codon_RF1):
  init_codon = initial_codon
  init_codon_stRNA = initial_codon_stRNA
  init_codon_RF1 = initial_codon_RF1
  init_cond = [init_codon,
              init_codon_stRNA,
              init_codon_RF1]
  time = list(range((3600))) #per minute for 24hrs
  return solve(init_cond, time)

def plot(time, to_plot):
  plt.plot(time, list([tup[0] for tup in to_plot]), label='codon')
  plt.plot(time, list([tup[1] for tup in to_plot]), label='codon:stRNA')
  plt.plot(time, list([tup[2] for tup in to_plot]), label='codon:RF1')
  plt.ylabel('Codons per State')
  plt.legend()
  plt.xlabel('Time (seconds)')
  plt.savefig('model.png')

results = []
params =[]
probabilities = []
counter = 0
for stRNA_bind in stRNA_bindings:
  for rf1_bind in rf1_bindings:
    for stRNA_unbind in stRNA_unbindings:
      for rf1_unbind in rf1_unbindings:
          for initial_codon in initial_codons:
            res = run(initial_codon, initial_codon_stRNA, initial_codon_RF1)
            results.append(res[1])
            params.append((stRNA_bind, rf1_bind, stRNA_unbind, rf1_unbind, initial_codon))
            probabilities.append(res[0])
            counter = counter + 1

res_prob = [[], [], [], [], []]
for i in range(len(probabilities)):
    if params[i][4] == 200.0:
        res_prob[0].append(probabilities[i])
    elif params[i][4] == 400.0:
        res_prob[1].append(probabilities[i])
    elif params[i][4] == 1000.0:
        res_prob[2].append(probabilities[i])
    elif params[i][4] == 2000.0:
        res_prob[3].append(probabilities[i])
    else:
        res_prob[4].append(probabilities[i])

maximums = [min(each) for each in res_prob]
minimums = [max(each) for each in res_prob]
avg = [sum(each)/len(each) for each in res_prob]

fig = plt.figure(figsize=(10, 10))

plt.fill_between([1, 2, 5, 10, 16], maximums, minimums)
plt.ylabel('Probability')
plt.yscale('log')
plt.legend()
plt.xlabel('Number of Stop Codons')
plt.savefig('prob_per_stop.png')

kanR_len = 819.
seconds_per_read = kanR_len/60.0
reads_per_hour = 3600.0/13.65
reads_per_hour_all_plasmids = reads_per_hour * 200

avg_fail_per_hour = [reads_per_hour_all_plasmids * each for each in avg]

fig = plt.figure(figsize=(10, 10))

plt.plot([1, 2, 5, 10, 16], avg_fail_per_hour)
plt.ylabel('Average Frequency of Failures per Hour')
plt.yscale('log')
plt.legend()
plt.xlabel('Number of Stop Codons')
plt.savefig('avg_fail_per_hour.png')
