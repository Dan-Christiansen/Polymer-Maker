import os
done_folders = [f for f in os.listdir(os.getcwd()) if os.path.isfile(f+'/nvt.gro')]
remaining_folders = [f for f in os.listdir(os.getcwd()) if f not in done_folders]
print("%6d Complete" % len(done_folders))
print("%6d Remaining" % len(remaining_folders))
print("%6d Total" % (len(done_folders)+len(remaining_folders)))
