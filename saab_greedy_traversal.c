typedef struct EdgeEls{
   int headNode;
   int tailNode;
   float score;
   EdgeEls *next;
} EdgeEl;

EdgeEls *listEdgeEl;

char *greedyGraphTraversal(float **graph_weights_matrix, int numNodes)
{
   char *tempResultSeq;
   float *alignScore;
   int overlapCount;
   float overlapScore;
   int **edgesPicked;
   char *resultSeq;

   float minScoreForAdding=15;
   int nodeChoosen;
   float maxValue=-1000;
   float minValueInList=-1000;
   int startNode;
   int stopNode;
   int prevStartNode=stopNode;	
   int prevStopNode=startNode;
   bool firstAdded, secondAdded;
   tempResultSeq=(char *) malloc(OEASeqMaxLen * sizeof(char));

   listEdgeEl=NULL;

   for (int count=0; count<numNodes; count++)
   {
      for (int count2=0; count2<numNodes; count2++)
      {
         if (graph_weights_matrix[count][count2]>maxValue && count!=count2)
         {
            addNewMaxEl(count, count2);
            maxValue=graph_weights_matrix[count][count2];

         }
         else if (graph_weights_matrix[count][count2] > maxValue - scoreVariance && count!=count2)
         {
            addNewEl(count, count2);
         }
      }
   }

   bestEdgeToPick(&startNode, &stopNode);

   read[startNode].mark=1;
   read[stopNode].mark=1;	
   alignScore=suffixPrefixAlignments(read[startNode].seq, read[stopNode].seq);
   overlapCount=-1;
   overlapScore=-1;
   for (int count=0; count<strlen(read[stopNode].seq); count++)
   {
      if (alignScore[count]>overlapScore)
      {
         overlapScore=alignScore[count];
         overlapCount=count;	
      }
   }

   free(alignScore);

   prevStartNode=startNode;
   prevStopNode=stopNode;



   resultSeq=addSeq(read[startNode].seq, read[stopNode].seq, overlapCount);

   freeList();
   listEdgeEl=NULL;

   do
   {
      maxValue=-1000;
      nodeChoosen=0;
      for (int count=0; count<numNodes; count++)
      {

         firstAdded=false;
         secondAdded=false;
         if (read[count].mark==0)
         {	      
            if (maxValue < graph_weights_matrix[count][prevStartNode])
            {  
               maxValue=graph_weights_matrix[count][prevStartNode];
               addNewMaxEl(count, prevStartNode);
               firstAdded=true;
            } 
            if (maxValue < graph_weights_matrix[prevStopNode][count])
            {	  
               maxValue=graph_weights_matrix[prevStopNode][count];
               addNewMaxEl(prevStopNode, count);
               secondAdded=true; 
            }
            if (firstAdded==false && maxValue-scoreVariance < graph_weights_matrix[count][prevStartNode])
            {

               addNewEl(count, prevStartNode);
            }
            if (secondAdded==false &&  maxValue - scoreVariance < graph_weights_matrix[prevStopNode][count])
            {

               addNewEl(prevStopNode, count);

            }
         }

      }
      bestEdgeToPick(&startNode, &stopNode);

      alignScore=suffixPrefixAlignments(read[startNode].seq, read[stopNode].seq);
      if (maxValue>minScoreForAdding)
      {
         if (startNode==prevStopNode)
         {
            free(alignScore);
            alignScore=suffixPrefixAlignments(read[prevStopNode].seq, read[stopNode].seq);
            overlapCount=-1;
            overlapScore=-100;
            for (int count=0; count<strlen(read[stopNode].seq)+1; count++)
            {
               if (alignScore[count]>overlapScore)
               {
                  overlapScore=alignScore[count];
                  overlapCount=count;
               }
            }
            strcpy(tempResultSeq, resultSeq);
            resultSeq=addSeq(resultSeq, read[stopNode].seq, overlapCount);

            if (strlen(resultSeq)>OEASeqMaxLen)
               return tempResultSeq;
            read[stopNode].mark=1;
            prevStopNode=stopNode;
            nodeChoosen=1;
            free(alignScore);

         }else if (stopNode==prevStartNode)
         {
            free(alignScore);
            alignScore=suffixPrefixAlignments(read[startNode].seq, read[prevStartNode].seq);
            overlapCount=-1;
            overlapScore=-1;
            for (int count=0; count<strlen(read[prevStartNode].seq)+1; count++)
            {
               if (alignScore[count]>overlapScore)
               {
                  overlapScore=alignScore[count];
                  overlapCount=count;
               }
            }
            strcpy(tempResultSeq, resultSeq);
            resultSeq=addSeq(read[startNode].seq,resultSeq, overlapCount);                                                       
            if (strlen(resultSeq)>OEASeqMaxLen)
               return tempResultSeq;
            read[startNode].mark=1;
            prevStartNode=startNode;
            nodeChoosen=1;
            free(alignScore);
         }
      }else free(alignScore);
      freeList();
      listEdgeEl=NULL;
   }while(maxValue>minScoreForAdding && strlen(resultSeq)<400);

   return(resultSeq);
}

char *addSeq(char *seq1, char *seq2, int overlapCount)
{
   char *resultSeq=(char *) malloc((strlen(seq1)+strlen(seq2)+1-overlapCount)*sizeof(char));
   for (int count=0; count<strlen(seq1); count++)
   {
      resultSeq[count]=seq1[count];
   }	
   for (int count=overlapCount; count<strlen(seq2); count++)
   {
      resultSeq[count+strlen(seq1)-overlapCount]=seq2[count];
   }
   resultSeq[(strlen(seq1)+strlen(seq2)-overlapCount)]='\0';
   return resultSeq;
}		

int addNewEl(int nodeStart, int nodeStop)
{
   EdgeEls *newEl;
   newEl= (EdgeEls *) malloc(sizeof(EdgeEls));
   newEl->headNode=nodeStart;
   newEl->tailNode=nodeStop;
   newEl->score=graph_weights_matrix[nodeStart][nodeStop];

   newEl->next=listEdgeEl;
   listEdgeEl=newEl;
}

int addNewMaxEl(int nodeStart, int nodeStop)
{
   EdgeEls *ptList;
   EdgeEls *headList;
   float maxValue;
   maxValue=graph_weights_matrix[nodeStart][nodeStop];  
   if (listEdgeEl!=NULL)
   {
      while (listEdgeEl!=NULL && listEdgeEl->score < maxValue-scoreVariance)
      {

         ptList=listEdgeEl->next;
         listEdgeEl=ptList;
      }

      while(listEdgeEl!=NULL && listEdgeEl->next!=NULL)
      {
         if (listEdgeEl->next->score < maxValue-scoreVariance)
         {
            ptList=listEdgeEl->next;
            listEdgeEl->next=ptList->next;
         }
         else 
            listEdgeEl=listEdgeEl->next;
      }
      addNewEl(nodeStart, nodeStop);

   }
   else
   {
      addNewEl(nodeStart, nodeStop);
   }
}


int bestEdgeToPick(int *source, int *dest)
{
   EdgeEls *ptrEdge;
   ptrEdge=listEdgeEl;
   float minValue=1000;
   int minValueSource=-1;
   int minValueDest=-1;
   while(ptrEdge!=NULL)
   {
      if (abs(read[ptrEdge->headNode].pos - read[ptrEdge->tailNode].pos)<minValue)
      {
         minValue=abs(read[ptrEdge->headNode].pos - read[ptrEdge->tailNode].pos);
         minValueSource=ptrEdge->headNode;
         minValueDest=ptrEdge->tailNode;
      }
      ptrEdge=ptrEdge->next;
   }
   *source=minValueSource;
   *dest=minValueDest;
}


void freeList(){
   EdgeEls *tmp;
   while (listEdgeEl != NULL){
      tmp = listEdgeEl->next;
      free(listEdgeEl);
      listEdgeEl = tmp;
   }
   listEdgeEl=NULL;
}
