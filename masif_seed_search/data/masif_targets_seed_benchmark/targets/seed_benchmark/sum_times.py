runtime = 0
with open('logs_site/all_times.txt') as f: 
    for line in f: 
        assert ('user' in line)
        usertime = line.split()[1]
        minutes = float(usertime.split('m')[0])
        seconds = float(usertime.split('m')[1].split('s')[0])
        runtime += minutes*60+seconds
print('runtime = {:.2f}'.format(runtime/60))
