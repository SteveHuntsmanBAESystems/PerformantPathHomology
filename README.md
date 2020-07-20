**README**

**OVERVIEW**
This is research code implementing "performant" path homology in MATLAB and written by M. Yutin. The calling file is pathhomology.m, which calls several other M-files. Note that this M-file has an identically named ancestor at https://github.com/SteveHuntsmanBAESystems/BasicPathHomology and upon which the present code drew from. The advantages of the ancestor code are its relative simplicity and permissive license; for serious applications, this code will perform better. It can benefit from MATLAB's symbolic toolbox to ensure accuracy that numerical linear algebra may not provide.

**CONTENTS**

    pathhomology.m              (wrapper; can use symbolic toolbox)
    pathhomology_allPaths.m     (computes allowed paths)
    pathhomology_bdryA.m        (computes unrestricted boundary operator)
    pathhomology_betti.m        (computes Betti numbers)
    pathhomology_homology.m     (computes homology representatives; can use symbolic toolbox)
    pathhomology_leafPaths.m    (reincorporates pruned leaves)
    pathhomology_lex.m          (generates lexicographical index)
    pathhomology_merge.m        (merges multiple homologies; can use symbolic toolbox)
    pathhomology_omega.m        (constructs the chain [path] complex; can use symbolic toolbox)
    pathhomology_parse.m        (preprocesses a digraph)
    README-PH.txt               (this file)
    removeloopsisos.m           (useful to preprocess digraphs; distributed under BSD-3)

**EXAMPLE** 
To run on an interesting example, try
        
        s = [1,1,2,2,2,3,4,4]; 
        t = [2,4,1,3,4,4,2,3]; 
        D = digraph(s,t); 
        ph = pathhomology(D,4)  % nontrivial homology in dimension 3

**TIPS**
Be careful about the size of input graphs and top dimension. Currently there are no provisions in the code to abort if things will get too big. This is left as an exercise.

Also, note that this code makes use of svds (via rank) if the symbolic mode is not chosen, and this sometimes leads to numerical errors, e.g., a negative Betti number. If the MATLAB symbolic toolbox is available, its use to confirm any findings (or fix any purportedly negative Betti numbers) is suggested.

**CITE**
If you use this code in your work, please cite it (naming only M. Yutin as author). We (Yutin and Huntsman) would also personally like to hear about your application, evaluation, etc.

**CONTACT**

    myutin@uscd.edu
    steve.huntsman@baesystems.com

**ACKNOWLEDGEMENTS**
This code is based upon work supported by BAE Systems FAST Labs. 

**LICENSE INFORMATION** 
"Performant" path homology files (including this one) are distributed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA. removeloopsisos.m is distributed under a BSD-3 license. All licenses are included in the files themselves, and this paragraph is for informational purposes only.

Copyright (c) 2020, BAE Systems. All rights reserved.
