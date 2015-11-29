function clusterMembershipProb = getClusterMembership(sampleCoords,embeddedPoints_train,trainData_genreID)
%{


%}


 totalGenres = length(unique(trainData_genreID));

 clusterMembershipMethod = 'kNN_of_kNN';

 switch clusterMembershipMethod
     
     case 'random'
         % % -----------------------------------------------------------------------
         % % Assigning Random membership probability Vector.
         probVec = rand(1,totalGenres);
         clusterMembershipProb = probVec/sum(probVec);
         
     case 'kNNvoting'
         
         % -----------------------------------------------------------------------
         % kNN - cluster-assignment by voting of kNN:
         dataSet = embeddedPoints_train;
         testSample = sampleCoords;
         
         [nbrs,dist_of_nbrs]=knnsearch(dataSet,testSample,'k',20,'distance','minkowski');
         
         genreIDs_nbrs = trainData_genreID(nbrs);
         probVec = zeros(1,totalGenres);
         
         for genreID = 1:totalGenres
             probVec(genreID) = sum(genreIDs_nbrs == genreID);
         end
         
   
     
     case 'kNN_of_kNN'
         % -----------------------------------------------------------------------
         % kNN - cluster-assignment by voting of kNN-of-kNN:
         
         dataSet = embeddedPoints_train;
         testSample = sampleCoords;
         
         [nbrs,dist_of_nbrs]=knnsearch(dataSet,testSample,'k',20,'distance','minkowski');
         
         nbr_coords = dataSet(nbrs,:);
         
         [nbrs_of_nbrs,dist_of_nbrs]=knnsearch(dataSet,nbr_coords,'k',20,'distance','minkowski');
         
         nbrs = reshape(nbrs,1,[]);
         nbrs_of_nbrs = reshape(nbrs_of_nbrs, 1, []);
         
         all_nbrs = [nbrs nbrs_of_nbrs];
         
         % -------------------------------------------------
         genreIDs_nbrs = trainData_genreID(nbrs);
         
         probVec = zeros(1,totalGenres);
         
         for genreID = 1:totalGenres
             probVec(genreID) = sum(genreIDs_nbrs == genreID);
         end

 end
 
 clusterMembershipProb = probVec / sum(probVec);

end