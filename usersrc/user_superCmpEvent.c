#include "user.h"

int user_superCmpEvent(superBurst *sbur,superCmpEvent *sevt)
{
    /* WARNING: do not alter things before this line */
    /*---------- Add user C code here ----------*/
    

    int cutWhichKilledEvent = SURVIVED;
#if ENABLE_EXCLUDE_BAD_BURSTS
    if( sbur->BadB.Dch!=0 || sbur->BadB.Phys!=0 )
    {
       // printf("47\n");
        cutWhichKilledEvent = 47;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif
    int numUntrackedClusters = 0;

    //we use this to store the indeces of the clusters that don't have a track
    int * tracklessClusters = NULL;

    user_lkrcalcor_SC (sbur,sevt,1); //lkr nonlinearity correction to all cluster energies

    /*** extract untracked clusters ****************************/
    for(int i = 0; (i < sevt->Ncluster);++i)
    {
        if(  sevt->cluster[i].iTrack == -1 && sevt->cluster[i].energy>3.0 )
            // iTrack = -1 if there is no associated track... greater than 3 GeV en to be a photon
        {
            ++numUntrackedClusters;
            tracklessClusters = (int *)realloc(tracklessClusters,
                                               numUntrackedClusters * sizeof(int));
            tracklessClusters[numUntrackedClusters - 1] = i;
        }

    }
#if ENABLE_DCH_LKR_VETO
    //DCH & Lkr veto:
    if( (numUntrackedClusters < 2) || ( sevt->Ntrack == 0) )
     {
      //  printf("72\n");
        cutWhichKilledEvent = 73;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif

    /***Min Bias Cut *************************/
    //this makes sure the data has the "minimum bias" flag set
    //it only rejects only the small fraction of the data which got through on the auto pass trigger
#if ENABLE_MIN_BIAS_CUT
    if( !((sevt->trigWord >> 11) & 1) )
    {
        cutWhichKilledEvent = 93;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif



    /*** LKr Quality cuts on trackless clusters **************************************************************/
    //clusters must be >2cm from dead cell, must have acceptable cluster status, and be in lkr acceptance

    int numGoodTracklessClusters=numUntrackedClusters;
    for(int i = 0;i<numUntrackedClusters;++i)
    {
        //the clusters positions at lkr face, corrected for projectivity:
        float clusterIPenetrationDepth = 20.8 + 4.3*logf(sevt->cluster[tracklessClusters[i]].energy);
        float lkrPlaneX = (sevt->cluster[tracklessClusters[i]].x + 0.136 + 0.00087*sevt->cluster[tracklessClusters[i]].y) *
                    (1+clusterIPenetrationDepth/10998);
        float lkrPlaneY = (sevt->cluster[tracklessClusters[i]].y + 0.300 - 0.00087*sevt->cluster[tracklessClusters[i]].x) *
                    (1+clusterIPenetrationDepth/10998);
        //printf("non corrected: %f, %f  corrected %f, %f \n",sevt->cluster[tracklessClusters[i]].x,sevt->cluster[tracklessClusters[i]].y,
          //          lkrPlaneX,lkrPlaneY);
        if( (sevt->cluster[tracklessClusters[i]].status > 4) || (sevt->cluster[tracklessClusters[i]].dDeadCell < 2) ||
                LKr_acc(sbur->nrun, lkrPlaneX, lkrPlaneY, 8)!=0   )
        {
            tracklessClusters[i]=-1;//if a cluster is bad- remove it by replacing its index with -1
            --numGoodTracklessClusters;
        }
    }
#if ENABLE_BASIC_QUALITY_CUTS
    if(numGoodTracklessClusters<2 && (cutWhichKilledEvent == SURVIVED))
    {
        cutWhichKilledEvent = 130;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif


    /*** Track Quality Cut **************************************************/
    // track must be in acceptance of DCHs
    //track must have less than 6% error on p measurement
    //track must be +ve - this only for p5 data which is +Kaons only.
    //must be within expected momentum range
    //must be within MUV acceptance
    //have a closest distance  of approach to run average beam position of less than 3.5 cm
    int numGoodTracks=sevt->Ntrack;//the number of potentially good tracks left
    int tracks[numGoodTracks], trackedClusters[numGoodTracks];//to store the index of good tracks, and the index of the associated clusters
    float tracksCDA[numGoodTracks],tracksTime[numGoodTracks], tracksChargedVertex[numGoodTracks][3];//to store the closest dist approach/vertex of each track
    float dzLkrDCH = Geom->Lkr.z - Geom->DCH.z;
    float dzMuv1DCH = Geom->Muv1.z - Geom->DCH.z;
    float dzMuv2DCH = Geom->Muv2.z - Geom->DCH.z;
    float dzMuv3DCH = Geom->Muv3.z - Geom->DCH.z;

    float beamPoint[3],beamVel[3];
    beamPoint[0]= abcog_params.pkxoffp;
    beamPoint[1] = abcog_params.pkyoffp;
    beamPoint[2] = 0.0;
    beamVel[0] = abcog_params.pkdxdzp;
    beamVel[1] = abcog_params.pkdydzp;
    beamVel[2] = 1.0;
    for(int i=0; i<sevt->Ntrack; i++)//check all tracks
    {
        float trackRadiusDCHb = sqrt(pow(sevt->track[i].bx,2)+pow(sevt->track[i].by,2));
        float trackRadiusDCH = sqrt(pow(sevt->track[i].x,2)+pow(sevt->track[i].y,2));
        float lkrPlaneX = sevt->track[i].x + dzLkrDCH*sevt->track[i].dxdz;
        float lkrPlaneY = sevt->track[i].y + dzLkrDCH*sevt->track[i].dydz;
        int trackCharge = sevt->track[i].q;
        float trackMomentum = sevt->track[i].p;
        float abCorrectedTrackMom = p_corr_ab(trackMomentum,trackCharge);//see http://goudzovs.web.cern.ch/goudzovs/ke2/selection.html

        float muv1x = sevt->track[i].x + dzMuv1DCH*sevt->track[i].dxdz;
        float muv1y = sevt->track[i].y + dzMuv1DCH*sevt->track[i].dydz;
        float muv2x = sevt->track[i].x + dzMuv2DCH*sevt->track[i].dxdz;
        float muv2y = sevt->track[i].y + dzMuv2DCH*sevt->track[i].dydz;
        float muv3x = sevt->track[i].x + dzMuv3DCH*sevt->track[i].dxdz;
        float muv3y = sevt->track[i].y + dzMuv3DCH*sevt->track[i].dydz;

        float chargedPartPoint[3],chargedPartVel[3];
        //for the charged particle here we work with the before magnetic field data so we can get
        //the vertex location
        chargedPartPoint[0] = sevt->track[i].bx;//x location of track in pre magnet DCH (DCHb)
        chargedPartPoint[1] = sevt->track[i].by;//y "
        chargedPartPoint[2] = Geom->DCH.bz;//z location of DCHb

        chargedPartVel[0] = sevt->track[i].bdxdz;//track moves in the dirc. (bdxdz,bdydz,1)
        chargedPartVel[1] = sevt->track[i].bdydz;
        chargedPartVel[2] = 1.0;

        float cda;//closest distance approach for each track
        float cdaVertex[3];
        //function from src/user.c:

        closap_(chargedPartPoint,beamPoint,chargedPartVel,beamVel,&cda,cdaVertex);

      //  fprintf(FP2,"%f\n",cda);
       // fprintf(FP1,"%f\n",cdaVertex[2]);

        tracksCDA[i] = cda;
        tracksChargedVertex[i][0]=cdaVertex[0];
        tracksChargedVertex[i][1]=cdaVertex[1];
        tracksChargedVertex[i][2]=cdaVertex[2];
        //get track time:
        if (sevt->track[i].hodstatus==2)
        {
                tracksTime[i] = sevt->track[i].hodTime;
        }
        else if (sevt->track[i].quality>0.9)
        {
                tracksTime[i]=sevt->track[i].time;
        }
        else if(sevt->track[i].hodstatus==1)
        {
                tracksTime[i] = sevt->track[i].hodTime;
        }
        else
        {
                tracksTime[i]=sevt->track[i].time;
        }
        //Quality cuts on the track:
        if(  (sevt->track[i].quality < 0.7) || LKr_acc(sbur->nrun, lkrPlaneX, lkrPlaneY, 8)!=0 || (trackRadiusDCH<14) || trackRadiusDCH>115 ||
               (trackRadiusDCHb<12) || trackRadiusDCHb>115 || ((sevt->track[i].perr/abCorrectedTrackMom)>0.06) || trackCharge!=1 || abCorrectedTrackMom<5 ||
                abCorrectedTrackMom>75 || muvAccept(muv1x,muv1y)!=1 || muvAccept(muv2x,muv2y)!=1 || muvAccept(muv3x,muv3y)!=1 || cda>5 ||
                (cdaVertex[2] < -2000) || (cdaVertex[2] > 9000) )
        {//we get rid of tracks that do not pass all of these
            --numGoodTracks;
            tracks[i]=-1;
            //printf("126\n");
        }
        else
        {//Cluster associated with track:
            int cluster = sevt->track[i].iClus;//this gives the index of the cluster or -1 if there is none

            //Cluster quality cuts:
            if( cluster>-1 && ( (sevt->cluster[cluster].status > 4) ||
                    (sevt->cluster[cluster].dDeadCell < 2) )  )
            {
                   --numGoodTracks;
                    tracks[i]=-1;
            }
            else
            {
                tracks[i]=i;//stores index of track
                trackedClusters[i]=sevt->track[i].iClus;//stores index of associated cluster
            }
        }
    }
#if ENABLE_BASIC_QUALITY_CUTS
    if( numGoodTracks<1 && (cutWhichKilledEvent == SURVIVED) )
    {
        cutWhichKilledEvent = 147;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif
 //   printf("tracks: %d tracks left: %d\n",sevt->Ntrack,numGoodTracks);


/*** Track and its Cluster Time Cut *******************************************/
//the timing of the track must be within 4ns of the associated lkr hit
#if ENABLE_TIMING_CUTS
    for(int i=0;i<sevt->Ntrack;++i)//for all tracks
    {
        if(tracks[i]>-1 && trackedClusters[i]>-1)//that haven't been discarded and have an associated cluster
        {//this wont work in MC (theres no timings)
             float dtTrackCluster =  tracksTime[i] - sevt->cluster[ trackedClusters[i] ].time;
             if (dtTrackCluster>4)
             {
                 tracks[i]=-1;
                 trackedClusters[i]=-1;
                 --numGoodTracks;
                 //printf("167\n");
             }
        }
    }

    if( numGoodTracks<1 && (cutWhichKilledEvent == SURVIVED) )
    {
        cutWhichKilledEvent = 278;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif




#if ENABLE_ENERGY_SHARING_CUTS
/****** Distance between in time gammas cut ****************************/
//gammas within 5ns must be more than 22cm apart to avoid energy sharing
    for(int i=0; i<numUntrackedClusters; i++)//for each gamma...
    {
        if(tracklessClusters[i] > -1)//that is good...
        {//using cluster projectivity corrections to get both positions at front face
            float clusterIPenetrationDepth = 20.8 + 4.3*logf(sevt->cluster[tracklessClusters[i]].energy);
            float clusterIx = (sevt->cluster[tracklessClusters[i]].x + 0.136 + 0.00087*sevt->cluster[tracklessClusters[i]].y) *
                    (1+clusterIPenetrationDepth/10998);
            float clusterIy = (sevt->cluster[tracklessClusters[i]].y + 0.300 - 0.00087*sevt->cluster[tracklessClusters[i]].x) *
                    (1+clusterIPenetrationDepth/10998);

            for(int s=0; s<numUntrackedClusters; s++)//comparing it and all the other gammas..
            {//that are within 5ns:
                float dtgigs = fabs(sevt->cluster[ tracklessClusters[i] ].time - sevt->cluster[ tracklessClusters[s] ].time);
                if((s!=i) && (dtgigs<5))//..and are also good... and not itself..
                {
                    float clusterSPenetrationDepth = 20.8 + 4.3*logf(sevt->cluster[tracklessClusters[s]].energy);
                    float clusterSx = (sevt->cluster[tracklessClusters[s]].x + 0.136 +
                                       0.00087*sevt->cluster[tracklessClusters[s]].y) * (1+clusterSPenetrationDepth/10998);
                    float clusterSy = (sevt->cluster[tracklessClusters[s]].y + 0.300 -
                                       0.00087*sevt->cluster[tracklessClusters[s]].x) * (1+clusterSPenetrationDepth/10998);

                    float distanceBetweenGammas_X = clusterSx - clusterIx;
                    float distanceBetweenGammas_Y = clusterSy - clusterIy;
                    float clusterClusterDist = sqrt(pow(distanceBetweenGammas_X,2) + pow(distanceBetweenGammas_Y,2));
                    if(clusterClusterDist<22)//if cluster i is within 22cm of another in time cluster...
                    {//then its a bad cluster
                        tracklessClusters[i]=-1;//if a cluster is bad- remove it by replacing its index with -1
                        --numGoodTracklessClusters;
                    }
                }
            }
        }
    }

/**** Cluster distance from all in time tracks cut **********************************/
//clusters must be > 22 cm away from any tracks that are within 10ns

    for(int n=0; n < sevt->Ntrack; ++n)//for all tracks...
    {
        for(int i=0; i < numUntrackedClusters; i++)//for all gammas
        {
            if((fabs( sevt->cluster[ tracklessClusters[i] ].time - tracksTime[n] ) < 10) )//track and cluster in time
            {
                float clusterIPenetrationDepth = 20.8 + 4.3*logf(sevt->cluster[tracklessClusters[i]].energy);
                float clusterIx = (sevt->cluster[tracklessClusters[i]].x + 0.136 + 0.00087*sevt->cluster[tracklessClusters[i]].y) *
                    (1+clusterIPenetrationDepth/10998);
                float clusterIy = (sevt->cluster[tracklessClusters[i]].y + 0.300 - 0.00087*sevt->cluster[tracklessClusters[i]].x) *
                    (1+clusterIPenetrationDepth/10998);

                float dzDchClusterZ = Geom->Lkr.z + clusterIPenetrationDepth - Geom->DCH.z;
                float TrackLkrPlaneX = sevt->track[n].x + dzDchClusterZ*sevt->track[n].dxdz;
                float TrackLkrPlaneY = sevt->track[n].y + dzDchClusterZ*sevt->track[n].dydz;

                float distanceFromTrack_X = clusterIx - TrackLkrPlaneX;
                float distanceFromTrack_Y = clusterIy - TrackLkrPlaneY;
                float clusterTrackDist = sqrt(pow(distanceFromTrack_X,2) + pow(distanceFromTrack_Y,2));

                if( (tracklessClusters[i] > -1) && clusterTrackDist<22)
                {
                    tracklessClusters[i]=-1;//if a cluster is bad- remove it by replacing its index with -1
                    --numGoodTracklessClusters;
                }
                if( tracks[n]>-1  && clusterTrackDist<22 )
                {
                    tracks[n]=-1;
                    --numGoodTracks;
                   // printf("266\n");
                }
            }
        }
    }

    if( ( (numGoodTracklessClusters<2) || (numGoodTracks<1) ) && (cutWhichKilledEvent == SURVIVED) )
    {
        cutWhichKilledEvent = 367;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif

/***Time difference between the tracks and untracked clusters ********************************************/
#if ENABLE_TIMING_CUTS
// the tracks be within 12 ns of atleast 2 photons
    for(int n=0; n<sevt->Ntrack; ++n)
    {
        if(tracks[n]>-1)
        {
            int numInTimePhotons=0;
            for(int i=0; i < numUntrackedClusters; i++)
            {
                if(tracklessClusters[i]>-1)
                {
                    float iClusterTime = sevt->cluster[ tracklessClusters[i] ].time;
                    if(  fabs(iClusterTime - tracksTime[n]) < 12 )//time cut 12ns
                    {
                        ++numInTimePhotons;
                    }
                }
            }
            if(numInTimePhotons<2)
            {
                tracks[n]=-1;
                --numGoodTracks;
               // printf("312\n");
            }
        }

    }
    //we require that the untracked cluster be within 12ns of at least one good track.
    for(int i=0; i < numUntrackedClusters; i++)
    {
        if( tracklessClusters[i] > -1)
        {
            int numInTimeTracks=0;
            float iClusterTime = sevt->cluster[ tracklessClusters[i] ].time;
            for(int n=0; n<sevt->Ntrack; ++n)
            {
                if(tracks[n] > -1)
                {
                    if(  (fabs(iClusterTime - tracksTime[n]) < 12) )//time cut 12ns
                    {
                        ++numInTimeTracks;
                    }
                }


            }
            if(numInTimeTracks<1)
            {
                tracklessClusters[i]=-1;
                --numGoodTracklessClusters;
            }
        }

    }

    if( ( (numGoodTracklessClusters<2) || (numGoodTracks<1) ) && (cutWhichKilledEvent == SURVIVED) )
    {
        cutWhichKilledEvent = 432;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif

 /****** time difference between the untrackedClusters *********************************/
 /*** selects pi0 gamma candidates ***/

//finds pairs of untracked clusters that arrive at lkr within 5ns of eachother, these are then candidate pi0 photons.
    float dtgigs;//time difference between gamma i and gamma s.
    int pi0GammaCandidatePairs[numGoodTracklessClusters*numGoodTracklessClusters][2];
    int numPi0GammaCandidatePairs=0;
    //this will probably be the best pi0 combination, however we will not assume so....
     //that is decided by z vertex location... if non are less than 2ns then we get rid of event
    for(int i=0; i<numUntrackedClusters; i++)//for each gamma...
    {
        if(tracklessClusters[i] > -1)//that is good...
        {
            for(int s=0; s<numUntrackedClusters; s++)
            //calcuate the time difference between it and all the other gammas..
            {
                if((s!=i) && (tracklessClusters[s] > -1) )//..that are also good... and not itself..
                {
                    dtgigs = fabs(sevt->cluster[ tracklessClusters[i] ].time -
                                        sevt->cluster[ tracklessClusters[s] ].time);

                    if(dtgigs<5)//if within 2ns of each other then we store them as pi0 gamma pair candidate
                    {

                        if(i<s)//makes sure we dont store same combination twice.
                        {

                            ++numPi0GammaCandidatePairs;
                            pi0GammaCandidatePairs[ numPi0GammaCandidatePairs-1 ][0] = tracklessClusters[i];
                            pi0GammaCandidatePairs[ numPi0GammaCandidatePairs-1 ][1] = tracklessClusters[s];

                            //now pi0Cands contains the cluster indices of the possible pi0 gamma pairs in each row
                        }
                    }

                }
            }
        }
    }
#if ENABLE_TIMING_CUTS

        if( (numPi0GammaCandidatePairs < 1) && (cutWhichKilledEvent == SURVIVED) )
            //if less than two - we dont have a pi0
        {
            cutWhichKilledEvent = 484;
            fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
#if BREAK_ON_FAILED_CUT
            return -1;
#endif
        }
#endif

/**** Track in time with Pi0 cut**********************************************/
// get rid of tracks that are not in time with one of the pi0 candidates.
  for(int n=0; n<sevt->Ntrack; ++n)
    {
        if(tracks[n]>-1)
        {
            int inTimeWithPi0=0;
            for(int i=0; i<numPi0GammaCandidatePairs; ++i)
            {
                float pi0time = (sevt->cluster[ pi0GammaCandidatePairs[i][0] ].time +
                                 sevt->cluster[ pi0GammaCandidatePairs[i][1] ].time )/2;
                if(   (fabs(pi0time - tracksTime[n] ) < 10) )
                {
                    ++inTimeWithPi0;
                }
            }
            if(inTimeWithPi0==0)
            {
                tracks[n]=-1;
                --numGoodTracks;
            }
        }

    }

/*** Find Blue Field Corrected charged Decay Vertex ***********************************************************/
//uses the algorith explained on http://goudzovs.web.cern.ch/goudzovs/ke2/selection.html
//cut on bf corrected zCharged and CDA
   // float tracksBfCorrChargedVertex[sevt->Ntrack][3];
    float zCharged;
    for(int i=0;i<sevt->Ntrack;++i)
    {
        if(tracks[i]>-1)
        {

            float chargedPartVel[3],chargedPartPoint[3],bfCorrChargedPartVel[3],bfCorrChargedPartPoint[3];;
            chargedPartVel[0] = sevt->track[tracks[i]].bdxdz;//track moves in the dirc. (bdxdz,bdydz,1)
            chargedPartVel[1] = sevt->track[tracks[i]].bdydz;
            chargedPartVel[2] = 1.0;
            chargedPartPoint[0] = sevt->track[tracks[i]].bx;//x location of track in pre magnet DCH (DCHb)
            chargedPartPoint[1] = sevt->track[tracks[i]].by;//y "
            chargedPartPoint[2] = Geom->DCH.bz;//z location of DCHb
            bfCorrChargedPartVel[0]=chargedPartVel[0];
            bfCorrChargedPartVel[1]=chargedPartVel[1];
            bfCorrChargedPartVel[2]=chargedPartVel[2];

            bfCorrChargedPartPoint[0]=chargedPartPoint[0];
            bfCorrChargedPartPoint[1]=chargedPartPoint[1];
            bfCorrChargedPartPoint[2]=chargedPartPoint[2];

            float decayVertex[3];
            decayVertex[0] = tracksChargedVertex[i][0];
            decayVertex[1] = tracksChargedVertex[i][1];
            decayVertex[2] = tracksChargedVertex[i][2];
            int trackCharge = 1;
            float trackMomentum = sevt->track[tracks[i]].p;
            float abCorrectedTrackMom = p_corr_ab(trackMomentum,trackCharge);
            blue_tack_(&trackCharge,&abCorrectedTrackMom,decayVertex,bfCorrChargedPartPoint,bfCorrChargedPartVel);
          //  printf("old slopes: %f ,%f new slopes: %f, %f \n",chargedPartVel[0],chargedPartVel[1],bfCorrChargedPartVel[0],bfCorrChargedPartVel[1]);
            //printf("old points: %f ,%f new points: %f, %f \n",chargedPartPoint[0],chargedPartPoint[1],bfCorrChargedPartPoint[0],bfCorrChargedPartPoint[1]);

            //float trackMidPointPosition[3];
            //trackMidPointPosition[2] = (Geom->DCH.bz - tracksChargedVertex[i][2])/2;
            //trackMidPointPosition[0] = chargedPartPoint[0] - trackMidPointPosition[2]*chargedPartVel[0];
            //trackMidPointPosition[1] = chargedPartPoint[1] - trackMidPointPosition[2]*chargedPartVel[1];
            //to get the coordinates of the decay vertex we use the charged tracks coordinates and velocity at DCHb (before magnet)
            //then extrapolate it back to ZneutralÍ
            float bfCorrVertex[3], bfCorrCda;
            closap_(bfCorrChargedPartPoint,beamPoint,bfCorrChargedPartVel,beamVel,&bfCorrCda,bfCorrVertex);
            //printf("cda: %f ---> %f\n",tracksCDA[i],bfCorrCda);
           // printf("zChar: %f ---> %f\n",tracksChargedVertex[i][2],bfCorrVertex[2]);

         //   fprintf(FP1,"%f\n",bfCorrCda);
        //fprintf(FP2,"%f\n",bfCorrVertex[2]);
            if( bfCorrCda>3 || bfCorrVertex[2]<-1600 || bfCorrVertex[2]>9000 )
            {
                tracks[i]=-1;
                --numGoodTracks;
                //tracksCDA[i]=bfCorrCda;
               // tracksBfCorrChargedVertex[i][0] = bfCorrVertex[0];
                //tracksBfCorrChargedVertex[i][1] = bfCorrVertex[1];
                //tracksBfCorrChargedVertex[i][2] = bfCorrVertex[2];
                zCharged=bfCorrVertex[2];

            }
        }
    }

#if ENABLE_Z_COORD_CUT
   if(numGoodTracks<1 && cutWhichKilledEvent == SURVIVED)
    {
        cutWhichKilledEvent = 572;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif


/*** Track Veto ************************************************************************/
//at this point the chance of having more than one good track is  small
//if it does happen we can just get rid of the event with little affect on the efficiency
    int iTrack;//=0;
    int iTrackedCluster;//=sevt->track[0].iClus;
   // float decayVertexClosestApproach[3];    //coordinates of decay
    //float iTrackCda;

   if(numGoodTracks!=1 && cutWhichKilledEvent == SURVIVED)
    {
        //printf("too many tracks\n");
        cutWhichKilledEvent = 534;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
    else
    {
        for(int i=0;i<sevt->Ntrack;++i)
        {
            if(tracks[i]>-1)
            {
                //fprintf(FP2,"%f\n",tracksChargedVertex[i][2]);
                //fprintf(FP1,"%f\n",tracksCDA[i]);
                iTrack = tracks[i];
                iTrackedCluster = sevt->track[i].iClus;

            }
        }

    }
  //  float zCharged = decayVertexClosestApproach[2];//Z coord of the vertex according to track
   // printf("%f %f\n",zCharged,iTrackCda);


/*** Here we calculate untracked cluster energies with CPD corrections *****************************/
//might need to make this do all good untracked clusters, depends if we decide to use that info or not.
    float tracklessClustersCorrectedEnergies[sevt->Ncluster];
    for(int i=0; i<numUntrackedClusters; ++i)
    {
        if(tracklessClusters[i]>-1)
        {
            // First find out to which cell is pointing the cluster hit (define CPDindex and CELLindex)
            // take lkr cluster position at front face
            double clusterPenetrationDepth = 20.8 + 4.3*logf(sevt->cluster[ tracklessClusters[i] ].energy);
            double clusterx = (sevt->cluster[ tracklessClusters[i] ].x + 0.136 + 0.00087*sevt->cluster[ tracklessClusters[i] ].y) *
                        (1+clusterPenetrationDepth/10998);
            double clustery = (sevt->cluster[ tracklessClusters[i] ].y + 0.300 - 0.00087*sevt->cluster[ tracklessClusters[i] ].x) *
                        (1+clusterPenetrationDepth/10998);
            int CELLindex;
            int CPDindex;
            //int * cpd_index=&CPDindex, *cell_index=&CELLindex;
            GetCpdCellIndex(clusterx, clustery, &CPDindex, &CELLindex);
            if( CELLindex==-1 || CPDindex==-1 )
            {
                printf("GetCpdCellIndex error\n");
                return -1;
            }
            // Ke3 E/p correction for each cell
            tracklessClustersCorrectedEnergies[ tracklessClusters[i] ] = sevt->cluster[  tracklessClusters[i]  ].energy / EopCorr[CPDindex][CELLindex];
           // printf("old: %f   new: %f \n", sevt->cluster[  tracklessClusters[i]  ].energy,tracklessClustersCorrectedEnergies[i]);
        }
    }

    /*** Calculates the decay vertex from cluster data ***********************************/
    /*** selects pi0 pair ***/

    float zDiffMin=1e308;
    int pi0pair=-1;
    float zPi,zNeutral;
    float pi0Photon1Energy,pi0Photon2Energy,pi0LkrEnergy;
    for(int i=0; i<numPi0GammaCandidatePairs; i++)
    {
        float gamma1Energy = tracklessClustersCorrectedEnergies[ pi0GammaCandidatePairs[i][0] ];
        float gamma2Energy = tracklessClustersCorrectedEnergies[ pi0GammaCandidatePairs[i][1] ];

        float gamma1PenetrationDepth = 20.8 + 4.3*logf(gamma1Energy);
        float gamma2PenetrationDepth = 20.8 + 4.3*logf(gamma2Energy);


        float gamma1LkrVertex[3],gamma2LkrVertex[3];

        gamma1LkrVertex[0] = (sevt->cluster[pi0GammaCandidatePairs[i][0]].x + 0.136 +
                           0.00087*sevt->cluster[pi0GammaCandidatePairs[i][0]].y) * (1+gamma1PenetrationDepth/10998);

        gamma1LkrVertex[1] = (sevt->cluster[pi0GammaCandidatePairs[i][0]].y + 0.300 -
                           0.00087*sevt->cluster[pi0GammaCandidatePairs[i][0]].x) * (1+gamma1PenetrationDepth/10998);

        gamma1LkrVertex[2] = Geom->Lkr.z;// + gamma1PenetrationDepth;

        gamma2LkrVertex[0] = (sevt->cluster[pi0GammaCandidatePairs[i][1]].x + 0.136 +
                           0.00087*sevt->cluster[pi0GammaCandidatePairs[i][1]].y) * (1+gamma2PenetrationDepth/10998);

        gamma2LkrVertex[1] = (sevt->cluster[pi0GammaCandidatePairs[i][1]].y + 0.300 -
                           0.00087*sevt->cluster[pi0GammaCandidatePairs[i][1]].x) * (1+gamma2PenetrationDepth/10998);

        gamma2LkrVertex[2] = Geom->Lkr.z;// + gamma2PenetrationDepth;

        float dispGamma1Gamma2[3];
        for(int j = 0; j < 3; ++j)
        {
            dispGamma1Gamma2[j] = gamma1LkrVertex[j] - gamma2LkrVertex[j];
        }

        float gamma1Gamma2Distance = f3vmag(dispGamma1Gamma2);
        //fprintf(FP1,"g1g2dist: %f gEn: %f %f  \n",gamma1Gamma2Distance,gamma1Energy,gamma2Energy);
        float dispPiLkr = gamma1Gamma2Distance*sqrt(gamma1Energy * gamma2Energy)/PI0_MASS;
        //disp of pi0 decay to lkr
        zPi = (Geom->Lkr.z) - dispPiLkr;//Z coord of vertex according to pi0
        //fprintf(FP2,"%f\n",zPi);
        //if( ((zPi > -1600) && (zPi < 9000)))
        {
            float zDiff = zCharged - zPi;
            if (zDiff<zDiffMin)
            {
                zDiffMin =zDiff;
                pi0pair = i;
                zNeutral = zPi;
                pi0Photon1Energy = gamma1Energy;
                pi0Photon2Energy = gamma2Energy;
                pi0LkrEnergy = pi0Photon1Energy +pi0Photon2Energy;
                

            }
        }
    }
  // fprintf(FP1,"%f\n",zNeutral);
#if ENABLE_Z_COORD_CUT
    if( ((zNeutral< -1600) || (zNeutral > 9000) || (pi0pair==-1) || pi0LkrEnergy<10 )&& (cutWhichKilledEvent == SURVIVED)  )
    {
        cutWhichKilledEvent = 727;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif
#if ENABLE_Z_COORD_CUT
    if( ( (zDiffMin< -8000) || (zDiffMin > 8000) ) && (cutWhichKilledEvent == SURVIVED)  )
    {
        cutWhichKilledEvent = 737;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);

#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif


/*** Decay Vertex calculation using z neutral and bf corrected slopes***/

    float chargedPartVel[3],chargedPartPoint[3],bfCorrChargedPartVel[3],bfCorrChargedPartPoint[3];;
    chargedPartVel[0] = sevt->track[iTrack].bdxdz;//track moves in the dirc. (bdxdz,bdydz,1)
    chargedPartVel[1] = sevt->track[iTrack].bdydz;
    chargedPartVel[2] = 1.0;
    chargedPartPoint[0] = sevt->track[iTrack].bx;//x location of track in pre magnet DCH (DCHb)
    chargedPartPoint[1] = sevt->track[iTrack].by;//y "
    chargedPartPoint[2] = Geom->DCH.bz;//z location of DCHb
    bfCorrChargedPartVel[0]=chargedPartVel[0];
    bfCorrChargedPartVel[1]=chargedPartVel[1];
    bfCorrChargedPartVel[2]=chargedPartVel[2];

    bfCorrChargedPartPoint[0]=chargedPartPoint[0];
    bfCorrChargedPartPoint[1]=chargedPartPoint[1];
    bfCorrChargedPartPoint[2]=chargedPartPoint[2];
    //using zn, first iteration vertex is
    float decayVertex[3];
    decayVertex[0] = chargedPartPoint[0]-chargedPartVel[0]*(chargedPartPoint[2]-zNeutral);
    decayVertex[1] = chargedPartPoint[1]-chargedPartVel[1]*(chargedPartPoint[2]-zNeutral);
    decayVertex[2] = zNeutral;
    int trackCharge = 1;
    float trackMomentum = sevt->track[iTrack].p;
    float abCorrectedTrackMom = p_corr_ab(trackMomentum,trackCharge);
    //get bf corrections;
    blue_tack_(&trackCharge,&abCorrectedTrackMom,decayVertex,bfCorrChargedPartPoint,bfCorrChargedPartVel);
    //second iteration vertex using bf corrections:
    decayVertex[0] = bfCorrChargedPartPoint[0]-bfCorrChargedPartVel[0]*(bfCorrChargedPartPoint[2]-zNeutral);
    decayVertex[1] = bfCorrChargedPartPoint[1]-bfCorrChargedPartVel[1]*(bfCorrChargedPartPoint[2]-zNeutral);
    decayVertex[2] = zNeutral;

    //now work out the displacement of decay vertex from beam:
    //first the beam position at Zn:
    float beamPositionZn[3];
    beamPositionZn[0] = beamPoint[0] + (zNeutral*beamVel[0]);//beamPoint/vel defined line 156
    beamPositionZn[1] = beamPoint[1] + (zNeutral*beamVel[1]);
    beamPositionZn[2] = zNeutral;//the z coord found from pi0 data

    float dispBeamDecayVertexX = beamPositionZn[0] - decayVertex[0];
    float dispBeamDecayVertexY = beamPositionZn[1] - decayVertex[1];

    float radialDistBeamDecayVertex= sqrt( pow(dispBeamDecayVertexY, 2) + pow(dispBeamDecayVertexX, 2) );

   // fprintf(FP1,"%f\n",dispBeamDecayVertexX);
    // fprintf(FP2,"%f %f\n",dispBeamDecayVertex[0],dispBeamDecayVertex[1]);
    //fprintf(FP2,"%f \n",radialDistBeamDecayVertex);

    //  fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
        //return 0;}

    //fprintf(FP2,"%f\n",zPi);
    //fprintf(FP1,"%f\n",zDiffMin);
    // printf("%f\n",zDiffMin);
#if ENABLE_Z_COORD_CUT
    //cut on radial decay distance
    if(radialDistBeamDecayVertex > 3 )
    {
        cutWhichKilledEvent = 795;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
#if BREAK_ON_FAILED_CUT
        return -1;
#endif
    }
#endif
//fprintf(FP1,"%f \n",radialDistBeamDecayVertex);
















/*** Kinematics Calulations *************************************************************************************************/
/***------------------------***/

    /*** Here we find the 2 pi0 gammas three momentum ***********************************************/
    //to work this out we need two points on the photons path
    //we use the points where each photon hits the Lkr:

    float gamma1Energy = pi0Photon1Energy;//to avoid rewriting code  below.
    float gamma2Energy = pi0Photon2Energy;

    float gamma1PenetrationDepth = 20.8 + 4.3*logf(gamma1Energy);
    float gamma2PenetrationDepth = 20.8 + 4.3*logf(gamma2Energy);
    float gamma1LkrVertex[3],gamma2LkrVertex[3];
    gamma1LkrVertex[0] = sevt->cluster[pi0GammaCandidatePairs[pi0pair][0]].x;
    gamma1LkrVertex[1] = sevt->cluster[pi0GammaCandidatePairs[pi0pair][0]].y;
    gamma1LkrVertex[2] = Geom->Lkr.z+gamma1PenetrationDepth;
    gamma2LkrVertex[0] = sevt->cluster[pi0GammaCandidatePairs[pi0pair][1]].x;
    gamma2LkrVertex[1] = sevt->cluster[pi0GammaCandidatePairs[pi0pair][1]].y;
    gamma2LkrVertex[2] = Geom->Lkr.z+gamma2PenetrationDepth;
    //and the second is the decay vertex...

    //we then have a vector of the photons paths after the decay
    float dispDecayVertexGamma1Vertex[3],dispDecayVertexGamma2Vertex[3];
    for(int i = 0; i < 3;++i)
    {
        dispDecayVertexGamma1Vertex[i] = gamma1LkrVertex[i] - decayVertex[i];
        dispDecayVertexGamma2Vertex[i] = gamma2LkrVertex[i] - decayVertex[i];
    }

    float distDecayVertexGamma1Vertex = f3vmag(dispDecayVertexGamma1Vertex);
    float distDecayVertexGamma2Vertex = f3vmag(dispDecayVertexGamma2Vertex);

    float gamma1VectMom[3],gamma2VectMom[3];
    for(int i = 0; i < 3; ++i)
    {
        //normalise then multiply by momentum magnitude (same as energy for a photon)
        gamma1VectMom[i] = gamma1Energy*(dispDecayVertexGamma1Vertex[i]/distDecayVertexGamma1Vertex);
        gamma2VectMom[i] = gamma2Energy*(dispDecayVertexGamma2Vertex[i]/distDecayVertexGamma2Vertex);
    }

    /****HERE WE RECONSTRUCT THE PI0 INVARIANT MASS************************/

    float pi0VectMom[3];
    for(int i = 0; i < 3; ++i)
    {
        pi0VectMom[i] = gamma1VectMom[i] + gamma2VectMom[i];
    }

    float pi0Energy = gamma1Energy + gamma2Energy;
    float pi0MomMag = f3vmag(pi0VectMom);

    float pi0FourMom[4]={ pi0Energy, pi0VectMom[0], pi0VectMom[1], pi0VectMom[2] };

    float pi0ReconstructedMass2 = f4vdot(pi0FourMom,pi0FourMom);
    float pi0ReconstructedMass = sqrt(pi0ReconstructedMass2);
    //sometimes reconstruction makes the mass imaginary
    //(this is bad so we ignore it and pretend it doesnt happen)
    //In c an imaginary number becomes "-nan" which has the property than -nan != -nan

    if( pi0ReconstructedMass != pi0ReconstructedMass)
    {
        printf("Error: pi0 mass reconstruction, line 593\n");
       // return -1;
    }

    //fprintf(FP1,"%f\n",pi0ReconstructedMass);



    /*** Here we calculate four momentum of additional (radiative) gammas *******************/
    //array to store the 4 mom of each additional photon
/*
    int numberOfAdditionalPhotons=0;
    if(numGoodTracklessClusters > 2)
    {
        float additionalGamma4Mom[numGoodTracklessClusters-2][4];

        for(int s=0; s<numGoodTracklessClusters; ++s)
        {//raditive gammas are those which did not get eliminated (set to -1) and are not the pi0
            if(tracklessClusters[s]>-1 && ( tracklessClusters[s]!=pi0GammaCandidatePairs[pi0pair][0] ||
                    tracklessClusters[s]!=pi0GammaCandidatePairs[pi0pair][1])   )
            {

                //now do the same as for the pi0 gammas
                float gammaLkrVertex[3];
                gammaLkrVertex[0] = sevt->cluster[s].x;
                gammaLkrVertex[1] = sevt->cluster[s].y;
                gammaLkrVertex[2] = Geom->Lkr.z;
                float dispDecayVertexGammaVertex[3];
                for(int i = 0; i < 3;++i)
                {
                     dispDecayVertexGammaVertex[i] = gammaLkrVertex[i] - decayVertex[i];
                }
                float distDecayVertexGammaVertex = f3vmag(dispDecayVertexGammaVertex);
                float gammaEnergy = tracklessClustersCorrectedEnergies[i];

                additionalGamma4Mom[numberOfAdditionalPhotons][0]=gammaEnergy;
                for(int i = 1; i < 4; ++i)
                {
                     additionalGamma4Mom[numberOfAdditionalPhotons][i] = gammaEnergy*(dispDecayVertexGammaVertex[i-1]/distDecayVertexGammaVertex);
                }
                ++numberOfAdditionalPhotons;
               // printf("radiative\n");
            }
        }
    }*/

/*    if(numGoodTracklessClusters > 2 && numberOfAdditionalPhotons==0 )
    {
        printf("------------------\n pi0 indices=%d, %d\n",pi0GammaCandidatePairs[pi0pair][0],pi0GammaCandidatePairs[pi0pair][1]);
        for(int i=0; i<numUntrackedClusters; ++i)
        {
            printf(" %d\t",tracklessClusters[i]);
        }
        printf("\n");
    }*/



/***Here we find 4 momentum of kaon(assumed to be beam average) and charged track **************************/
    //first we find momentum of kaon, assumed to be beam average.
    //this is really dodgy because we've already shown that abcog_params lies to us!
    //need to find out how they do it in "Reboot" slides.
    float beamVelocity[3];
    beamVelocity[0] = abcog_params.pkdxdzp;
    beamVelocity[1] = abcog_params.pkdydzp;
    beamVelocity[2] = 1.0;
    float beamVelMag = f3vmag(beamVelocity);

    float beamMomMag =  p_corr_ab(abcog_params.pkp,1);//alpha-Beta correction
    //E = (p^2 + m^2)^1/2
    float beamKaonEnergy = sqrt( pow(beamMomMag,2) + pow(abcog_params.mkp,2) );
    float kaonEnergy = beamKaonEnergy;//sqrt(f3vmag2(beamVectMom) + _POWER2(abcog_params.mkp) );

    //track:
    float chargedPartMomMag = p_corr_ab(sevt->track[iTrack].p,1);
    //using chargedPartVel declared line 592
    float chargedPartVelMag = f3vmag(chargedPartVel);



/*** He we find track energy with CPD Corrections ****************************************************/
    float chargedPartEnergy;
    if ( iTrackedCluster > -1)
    {
        // First find out to which cell is pointing the deflected track (define CPDindex and CELLindex)
        //find track coords at lkr face:
        float dzDchClusterZ = Geom->Lkr.z - Geom->DCH.z;
        double trkatlkr[2];
        trkatlkr[0] = sevt->track[iTrack].x + dzDchClusterZ*sevt->track[iTrack].dxdz;
        trkatlkr[1] = sevt->track[iTrack].y + dzDchClusterZ*sevt->track[iTrack].dydz;

        int CELLindex;
        int CPDindex;

        GetCpdCellIndex(trkatlkr[0] , trkatlkr[1], &CPDindex, &CELLindex);

        // Now that you know the cell hit by the track, correct the energy for the track cluster (is there is one associated to the track)
        chargedPartEnergy = sevt->cluster[iTrackedCluster].energy / EopCorr[CPDindex][CELLindex];  // Ke3 E/p correction for each cell
    }
    else
    {
        chargedPartEnergy = 0;
    }

    float chargedPartEPRatio = chargedPartEnergy / chargedPartMomMag;

    fprintf(FP1,"%f\n",chargedPartEPRatio);

    float chargedPartVectMom[3], beamVectMom[3];

    for(int i = 0; i < 3; ++i)
    {
        chargedPartVectMom[i] = chargedPartMomMag*chargedPartVel[i]/chargedPartVelMag;
        beamVectMom[i] = beamMomMag * beamVelocity[i] / beamVelMag;
    }
    float chargedPart4Mom[4] = {chargedPartEnergy, chargedPartVectMom[0], chargedPartVectMom[1], chargedPartVectMom[2] };
    float beam4Mom[4] = { beamKaonEnergy, beamVectMom[0], beamVectMom[1], beamVectMom[2]  };


    /*****Missing (three) Momentum Calculation ******************************************/

    float detectedMom[3];
    for(int i = 0; i < 3; ++i)
    {
        detectedMom[i] = pi0VectMom[i] + chargedPartVectMom[i];
    }

    float missingMom[3];
    for( int i = 0; i < 3; ++i)
    {
        missingMom[i] = beamVectMom[i] - detectedMom[i];
    }

    float missingMomMag = f3vmag(missingMom);

    /****HERE WE ASSUME THE CHARGED PARTICLE IS A PI+ AND CALCULATE THE KAON MASS*************/
    // see arXiv:hep-ex/0702015v2 §4.2
    //see the note in user.h about naming the pi+ Pi
    float pi0EnergyNonLkrMeasure = sqrt( _POWER2(PI0_MASS) + _POWER2(pi0MomMag));
    float pi1Energy = sqrt(_POWER2(PI1_MASS) + _POWER2(chargedPartMomMag));
    float pi0Pi1Mass = sqrt(_POWER2(pi0EnergyNonLkrMeasure + pi1Energy) - f3vmag2(detectedMom));
    /*** same, assuming electron, and muon**************************************************/
    float eEnergy = sqrt(_POWER2(ELECTRON_MASS) + _POWER2(chargedPartMomMag));
    float pi0ElectronMass =  sqrt(_POWER2(pi0EnergyNonLkrMeasure + eEnergy) - f3vmag2(detectedMom));
    float muEnergy = sqrt(_POWER2(MUON_MASS) + _POWER2(chargedPartMomMag));
    float pi0MuonMass =  sqrt(_POWER2(pi0EnergyNonLkrMeasure + muEnergy) - f3vmag2(detectedMom));

   // fprintf(FP2,"%f\n",pi0Pi1Mass);


    /***find assoiciated muon if present*************************************/

    //This is -1 if there is no muon, apparently with high reliability
    //the muon structure is filled by the murec0902 routine
    // source code: /afs/cern.ch/user/g/goudzovs/offline/compact/compact-7.3/compact/rlib/anasrc/murec0902.c
    //it leaves you to check which planes the was a hit in, encoded in the status variable, we require hits in all 3 or just 1 and 2.
    //it leaves you to test the time resultion, we require within 3.5 ns
    int iMuon = sevt->track[iTrack].iMuon;//index of associated muon
    int muon;//if we
    if(iMuon == -1)
    {
        muon = 0;
    }
    else
    {
        float timeBetweenTrackAndMuon = fabs( tracksTime[iTrack] - sevt->muon[iMuon].time);
        if( (sevt->muon[iMuon].status==1 || sevt->muon[iMuon].status==2) && timeBetweenTrackAndMuon<=3.5 )
        {
            fprintf(FP2,"%f\n", chargedPartEPRatio);
            muon = 1;//then we have a muon
        }
    }

    /****Here we calculate the transverse momenta ********************************/
    float beamDirection[3];
    for(int i = 0;i<3;++i)
    {
        beamDirection[i] = beamVelocity[i]/beamVelMag;
    }
    //for the pi0:
    float pi0BeamDirMomMag = f3vdot(pi0VectMom,beamDirection);
    float pi0TransMom[3];
    for(int i = 0;i<3;++i)
    {
        pi0TransMom[i] = pi0VectMom[i] - pi0BeamDirMomMag*beamDirection[i];
    }
    float pi0TransMomMag = f3vmag(pi0TransMom);
   // fprintf(FP1,"%f\n",pi0TransMomMag);
    //the track:
    float trackBeamDirMomMag = f3vdot(chargedPartVectMom,beamDirection);
    float trackTransMom[3];
    for(int i = 0;i<3;++i)
    {
        trackTransMom[i] = chargedPartVectMom[i] - trackBeamDirMomMag*beamDirection[i];
    }
    float trackTransMomMag = f3vmag(trackTransMom);

    //event:
    float eventTransMom[3];
    for(int i = 0;i<3;++i)
    {
        eventTransMom[i] = trackTransMom[i] + pi0TransMom[i];
    }
    float eventTransMomMag = f3vmag(eventTransMom);
   // fprintf(FP2,"%f\n",eventTransMomMag);
       
    ++nEvents;

    /**** Here we I.D. The Track ********************************************************/
    int trackID;
    float missingMassSquared;
    //positron selection criteria
    if( (pi0Pi1Mass<0.4772 || pi0Pi1Mass>0.5102) && trackTransMomMag<0.2 && chargedPartEPRatio>0.9 && muon==0 && chargedPartMomMag<35 &&
            pi0ElectronMass<0.425  )
    {
        trackID = ELECTRON;
        missingMassSquared = pow(kaonEnergy,2) + pow(pi0Energy,2) + pow(ELECTRON_MASS, 2) + pow(chargedPartMomMag,2) +
            2*pi0Energy*sqrt( pow(ELECTRON_MASS, 2) + pow(chargedPartMomMag,2) ) - 2*kaonEnergy*pi0Energy -
            2*kaonEnergy*sqrt( pow(ELECTRON_MASS, 2) + pow(chargedPartMomMag,2) ) -
            pow(missingMomMag,2);
        fprintf(ke3FP,"%f\n",missingMassSquared);
        if(missingMassSquared>-0.0120 && missingMassSquared <0.012)
            ++ke3Count;
    }
    //muon selection criteria
    else if( (pi0Pi1Mass<0.4772 || pi0Pi1Mass>0.5102) && trackTransMomMag<0.2 && muon==1 && pi0Energy<40 && pi0MuonMass< 0.38 && 
            chargedPartMomMag>10 && chargedPartMomMag<40)
    {
        trackID = MUON;
        missingMassSquared = pow(kaonEnergy,2) + pow(pi0Energy,2) + pow(MUON_MASS, 2) + pow(chargedPartMomMag,2) +
            2*pi0Energy*sqrt( pow(MUON_MASS, 2) + pow(chargedPartMomMag,2) ) - 2*kaonEnergy*pi0Energy -
            2*kaonEnergy*sqrt( pow(MUON_MASS, 2) + pow(chargedPartMomMag,2) ) -
            pow(missingMomMag,2);

        fprintf(km3FP,"%f\n",missingMassSquared);
        if(missingMassSquared>-0.01 && missingMassSquared <0.01)
            ++km3Count;
    }
    //pi+ selection criteria
    else if( (pi0Pi1Mass>0.4772 && pi0Pi1Mass<0.5102) && trackTransMomMag<0.215 && chargedPartEPRatio<0.9 && chargedPartMomMag>10 && 
            chargedPartMomMag<50)
    {
        trackID = PIPLUS;
        missingMassSquared = pow(kaonEnergy,2) + pow(pi0Energy,2) + pow(PI1_MASS, 2) + pow(chargedPartMomMag,2) +
            2*pi0Energy*sqrt( pow(PI1_MASS, 2) + pow(chargedPartMomMag,2) ) - 2*kaonEnergy*pi0Energy -
            2*kaonEnergy*sqrt( pow(PI1_MASS, 2) + pow(chargedPartMomMag,2) ) -
            pow(missingMomMag,2);

        fprintf(k2piFP,"%f\n",missingMassSquared);
        if(missingMassSquared>-0.0025 && missingMassSquared <0.001)
            ++k2piCount;
    }
    else
    {
        trackID = UNIDENTIFIED;
        ++numUnidentified;
                cutWhichKilledEvent = 1008;
        fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);
        return -1;
    }



    /****HERE WE DO THE MISSING MASS CALCULATION***********************************/

    //float missingEnergy = kaonEnergy - pi0Energy - chargedPartEnergy;
    //we worked out the missing momentum above

  //  float missingMass2 = _POWER2(missingEnergy) - _POWER2(missingMomMag);

    //new method avoiding using Etrack since this is not accurate for pi+ and muon tracks







    //fprintf(FP2,"%f\n",missingMass2);


#if ENABLE_OUTPUT

   // fprintf(realFP,"%f %f %f %f %f %d %f\n",pi0ReconstructedMass, pi0Pi1Mass, missingMomMag, chargedPartEPRatio,missingMass2,muon,pi0TransMomMag);
    //  fprintf(realFP,"%f %f %f\n",decayVertex[0],decayVertex[1],decayVertex[2]);
#endif

    free(tracklessClusters);

    fprintf(cutWhichKilledEventFP,"%d\n",cutWhichKilledEvent);


   // nuserevt++;

    /*----------- End of user C code -----------*/
    return 0;
}
