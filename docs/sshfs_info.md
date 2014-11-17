# How to mount the FAT filesystem:
```bash
[jah5 ~/Desktop/Classes/COS597D/cos597project]$ ssh -f med_impute -L 2222:fat:22 -N
[jah5 ~/Desktop/Classes/COS597D/cos597project]$ sshfs -p 2222 jah5@localhost:/memex/jah5/notebooks ../fat_mount/
```
