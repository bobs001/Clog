Git Workflow Summary
Syncing Central Repo with Local Repo

Setting It Up (only do this the initial time)

    Find & copy Central Repo URL
    git remote add upstream https://github.com/bobs001/Clog.git

After Initial Set Up

    Update your Local Repo & Push Changes
        git pull upstream master - pull down any changes and sync the local repo with the central repo
        make changes, git add and git commit
        git push origin master - push your changes up to your fork
        Repeat

- 3/3/2020
  Cloned Clog repo
  $ git clone https://github.com/bobs001/Clog.git Clog

- 3/3/2020
  I want to push changes from local Clog back to https://github.com/bobs001/Clog
  $ cd Clog
  $ git remote add upstream https://github.com/bobs001/Clog.git [did not quite work]
  Edited .git/config by hand.  
  $ git push


