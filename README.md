# PLasserreZuber_labBook  

Internal bioinfo projects as support to research teams of the INRAe UMR GDEC.  
From /home/palasser/projects/soutien_bioinfo on HPC2 cluster  


### gitlab    

```bash
## local git repository
cd /home/palasser/projects/soutien_bioinfo
git init
git add */bin/*
git commit -m "initial commit"
git status
## remote git repository
git remote add sb https://forgemia.inra.fr/gdec-bioinfo/plasserrezuber_labbook.git
git remote -v
git push -u sb --all
## merge request
## local git repository
git branch -M main
git branch -a
git fetch sb main
git pull sb main
```


### Support  
pauline.lasserre-zuber@inrae.fr
