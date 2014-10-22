## For running ipython notebook on the fat server and editing on your laptop

(1) connect via `ssh` to fat via : `ssh med_impute; ssh fat`
(1a) do `source ~/.bashrc`
(2) Start ipython notebook in the `notebooks` directory: `cd notebooks; ipython notebook --profile=jahserver`
(3) On your laptop, in another shell, forward port 1234 on FAT to your laptop: `ssh -f -L 1234:fat:1234 jah5@mmx.cs.princeton.edu -N`
(4) Open a browser and go to `localhost:1234`
