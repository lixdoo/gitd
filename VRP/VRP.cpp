// VRP.cpp : �������̨Ӧ�ó������ڵ㡣
//

//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <time.h>
#include <math.h>

struct Point
{
	double x; 
	double y;
};

const int NUM_ANT=60;
const int NUM_RUN=5;
const int ALPHA=2;
const int BETA=2;
const int GAMA=2;
const double ROU=0.5;
const double RANGE=20;
const int Q1=20;
const int Q2=20;
const int Q3=20;
const double SMALL=0.02;
const double VERY_SMALL=0.00001;
const int MAX=0x7fffffff;

double Random(double from, double to);
double FindDistance(Point a, Point b);
double FindPower(double in, int power);
int FindNext(double prob, int nClient, double* aProb/*, bool* aCntn*/);
int FindNextAnt(int from, int& nFind, const int** aTaboo, bool* aSuccess, 
				bool* aFound, int* aSucAnt, const int nClient, int& nAnt);
void FindSolution(int from, int nFind, const int** aTaboo, bool* aSuccess, double** aDeltaTaoP, int& nPath,
				  double& opt, int* optAnt, int& nOptAnt, int nAnt, double** aDstn, bool* aFound, int* aSucAnt, const int nClient);




int main()
{
	FILE* pIn=NULL;
	FILE* pResult=NULL;
	int i, j, k, nRun, nTime, before;			//, nLocalRun;
	int nClient, load, next, nPath, nOptAnt;
	Point* aPos=NULL;
	double** aDstn=NULL;			// ����
	double** aDeltaDstn=NULL;
	int* aDmnd=NULL;				// �ͻ�����
	int minDmnd=MAX;
	bool* aCntn=NULL;				// �Ƿ��ѻص���������
	int** aTaboo=NULL;			
	int* aWeight=NULL;				// ��ǰ����
	int nEnd;								// ���ε�����������������������
//	int nZero;								// �ڵ���������Ϊ0�Ľڵ���
	double** aTao=NULL;
	double** aDeltaTaoP=NULL;
	double deltaTao;
	double* aNAnt=NULL;					// �õ�������  �ø���������ֹ���
	double** aNTrans=NULL;				// ��iת�Ƶ�j��������  �ø���������ֹ���
	int* aLastPath=NULL;			// ��¼��һ���ߵ��ĵ�
	double* aProb=NULL;
	bool* aSuccess=NULL;			// �������ߵ�·���Ƿ��ǿ���·��
	bool* aFound=NULL;
	int* aSucAnt=NULL;
	int* optAnt=NULL;
	double prob;
	double sumProb;
	double opt, lastopt=(double)MAX;
	int mode, nLine;

	srand((unsigned)time(NULL));
	
	pIn=fopen("D:\\DataIn.txt", "r");
	if (!pIn)
	{
		return -1;
	}
	fscanf(pIn, "%d\n%d\n", &mode, &nClient);

	aPos=new Point [nClient+1];
	aPos[0].x=0.0;
	aPos[0].y=0.0;
	aDstn=new double* [nClient+1];
	aDeltaDstn=new double* [nClient+1];
	aTao=new double* [nClient+1];
	aDeltaTaoP=new double* [nClient+1];
//	aDeltaTao=new double* [NUM_ANT];
	aNAnt=new double [nClient+1];
	aNTrans=new double* [nClient+1];
	aProb=new double [nClient+1];
	aFound=new bool [nClient+1];
	for (i=0; i<=nClient; i++)
	{
		aDstn[i]=new double [nClient+1];
		aDeltaDstn[i]=new double [nClient+1];
		aTao[i]=new double [nClient+1];
		aDeltaTaoP[i]=new double [nClient+1];
//		aDeltaTao[i]=new double [nClient+1];
		aNTrans[i]=new double [nClient+1];
	}
	aDmnd=new int [nClient+1];
	aDmnd[0]=0;

	if (mode==0)				// ��Ҫ���ͻ�����λ�õĳ�ʼ������
	{
		for (i=1; i<=nClient; i++)
		{
			aPos[i].x=Random(-1*RANGE, RANGE);
			aPos[i].y=Random(-1*RANGE, RANGE);
		}
	}
	else
	{
		for (i=1; i<=nClient; i++)
		{
			fscanf(pIn, "%lf\t%lf\n", &aPos[i].x, &aPos[i].y);
		}
	}
	for (i=0; i<=nClient; i++)
	{
		aDstn[i][i]=0.0;
		for (j=i+1; j<=nClient; j++)
		{
			aDstn[i][j]=FindDistance(aPos[i], aPos[j]);
			aDstn[j][i]=aDstn[i][j];
		}
	}

	for (i=0; i<=nClient; i++)
	{
		for (j=i; j<=nClient; j++)
		{
			aDeltaDstn[i][j]=FindPower(aDstn[i][0]+aDstn[0][j]-aDstn[i][j]+SMALL, GAMA);
			aDeltaDstn[j][i]=aDeltaDstn[i][j];
		}
	}
	
	fscanf(pIn, "%d", &load);
	for (i=1; i<=nClient; i++)
	{
		fscanf(pIn, "%d ", aDmnd+i);
		if (minDmnd>aDmnd[i])
		{
			minDmnd=aDmnd[i];
		}
	}

	fclose(pIn);

	aTaboo=new int* [NUM_ANT];
	aCntn=new bool [NUM_ANT];
	aLastPath=new int [NUM_ANT];
	aSuccess=new bool [NUM_ANT];
	optAnt=new int [NUM_ANT];
	aSucAnt=new int [NUM_ANT];
	for (i=0; i<NUM_ANT; i++)
	{
		aTaboo[i]=new int [nClient+1];
	}
	aWeight=new int [NUM_ANT];
	

	for (i=0; i<=nClient; i++)
	{
//		aNAnt[i]=0.0;
		for (j=0; j<=nClient; j++)
		{
			aTao[i][j]=SMALL;			// ��Ϊ����һ������ֵ
			//				aDeltaTao[i][j]=0.0;
//			aNTrans[i][j]=0.0;
		}
	}
	// ��ʼ��Ⱥ�㷨
	for (nRun=0; nRun<NUM_RUN; nRun++)
	{
		// ��ʼ��

		nPath=0;
		opt=(double)MAX;

		for (i=0; i<=nClient; i++)
		{
			aNAnt[i]=0.0;
			for (j=0; j<=nClient; j++)
			{
//				aTao[i][j]=1.0;			// ��Ϊ����һ������ֵ
//				aDeltaTao[i][j]=0.0;
				aNTrans[i][j]=0.0;
			}
		}
		aNAnt[0]=NUM_ANT;
		for (i=0; i<NUM_ANT; i++)
		{
			aWeight[i]=load;
			aCntn[i]=true;
			for (j=0; j<=nClient; j++)
			{
				aTaboo[i][j]=-1;
			}
		}

		// �ߵ�һ��
		for (i=0; i<NUM_ANT; i++)
		{
			/*
			sumProb=0.0;
			for (j=1; j<=nClient; j++)
			{
				sumProb+=FindPower(aTao[0][j], ALPHA)*FindPower(1.0/aDstn[0][j], BETA)*aDeltaDstn[0][j];
			}
			aProb[0]=0.0;
			for (j=1; j<=nClient; j++)
			{
				aProb[j]=aProb[j-1]+FindPower(aTao[0][j], ALPHA)*FindPower(1.0/aDstn[0][j], BETA)
					*aDeltaDstn[0][j]/sumProb;
			}
	
			prob=Random(0.0, 1.0);
			*/

			next=(int)Random(1.0, (double)(nClient+1.0));	//FindNext(prob, nClient, aProb, aCntn);
			aLastPath[i]=next;

			aTaboo[i][next]=0;
			aWeight[i]-=aDmnd[next];
			aNAnt[next]+=1.0;
			aNTrans[0][next]+=1.0;
		}
		// ������Ϣ��
		
		for (i=1; i<=nClient; i++)
		{
//			deltaTao=0.0;
//			for (j=0; j<NUM_ANT; j++)
//			{
//				if (aTaboo[j][i]==0)
//				{
//					deltaTao+=Q1*(1.0-aNTrans[0][i]/aNAnt[0]);
//				}
//			}
			deltaTao=Q1*(1.0-aNTrans[0][i]/aNAnt[0]);
			aTao[0][i]=aTao[0][i]*(1-ROU)+deltaTao;
		}
		
		// �������ɲ�
		nEnd=0;
		nTime=1;
//		nZero=0;
		aNAnt[0]=0;
		for (i=0; i<NUM_ANT; i++)
		{
			aSuccess[i]=true;
		}

		while (nEnd<NUM_ANT)	// && nZero<nClient
		{
			for (i=0; i<NUM_ANT; i++)
			{
				if (!aCntn[i])
				{
					continue;
				}

				sumProb=0.0;
				for (j=0; j<=nClient; j++)
				{
					if (aWeight[i]-aDmnd[j]>=0 && aTaboo[i][j]<0)
					{
						sumProb+=FindPower(aTao[aLastPath[i]][j], ALPHA)
							*FindPower(1.0/(aDstn[aLastPath[i]][j]+SMALL), BETA)*aDeltaDstn[aLastPath[i]][j];
					}
				}
				//aProb[0]=0.0;
				aProb[0]=FindPower(aTao[aLastPath[i]][0], ALPHA)*
					FindPower(1.0/aDstn[aLastPath[i]][0], BETA)*aDeltaDstn[aLastPath[i]][0]/(sumProb+SMALL);
				if (aProb[0]<VERY_SMALL)
				{
					aProb[0]=0.0;
				}
				for (j=1; j<=nClient; j++)
				{
					if (aWeight[i]-aDmnd[j]>=0 && aTaboo[i][j]<0)
					{
						aProb[j]=aProb[j-1]+FindPower(aTao[aLastPath[i]][j], ALPHA)
							*FindPower(1.0/(aDstn[aLastPath[i]][j]+SMALL), BETA)*aDeltaDstn[aLastPath[i]][j]/(sumProb+SMALL);
					}					
					else
					{
						aProb[j]=aProb[j-1];
					}
				}

				prob=Random(-0.6, aProb[nClient]);
				next=FindNext(prob, nClient, aProb/*, aCntn*/);

//				aNAnt[aLastPath[i]]--;
				/*
				if (aWeight[i]-aDmnd[next]<0)
				{
					if (aWeight[i]-minDmnd<0)
					{
						aCntn[i]=false;
						nEnd++;
						aNTrans[aLastPath[i]][0]++;
					}
					continue;
				}
				*/
				aWeight[i]-=aDmnd[next];

				aNAnt[next]+=1.0;
				aNTrans[aLastPath[i]][next]+=1.0;
//				aTaboo[i][next]=nTime;

				if (next==0)
				{
					aCntn[i]=false;
					nEnd++;
					if (aTaboo[i][0]<0)
					{
						aTaboo[i][0]=nTime;
					}
					continue;
				}
//				if (aTaboo[i][next]>-1)
//				{
//					aTaboo[i][0]=nTime;           // ??????????????
//					aSuccess[i]=false;               // ??????????????
//					nEnd++;
//				}
				else if (aTaboo[i][0]<0)
				{
					aTaboo[i][next]=nTime;
				}
				
				aLastPath[i]=next;
			}

			for (k=0; k<=nClient; k++)
			{
				for (i=0; i<=nClient; i++)
				{
					deltaTao=0.0;
//					for (j=0; j<NUM_ANT; j++)
//					{
//						if (aTaboo[j][i]==nRun)
//						{
//							deltaTao+=Q1*(1.0-aNTrans[k][i]/aNAnt[k]);
//						}
//					}
					deltaTao=Q1*(1.0-aNTrans[k][i]/aNAnt[k]);
					aTao[k][i]=aTao[k][i]*(1-ROU)+deltaTao;
				}
			}

			nTime++;
		}

		// �ҿ��н�
		for (i=0; i<=nClient; i++)
		{
			aFound[i]=false;
			for (j=0; j<=nClient; j++)
			{
				aDeltaTaoP[i][j]=0.0;
			}
		}
		
		FindSolution(0, 0, (const int**)aTaboo, aSuccess, aDeltaTaoP, nPath, opt, 
			optAnt, nOptAnt, 0, aDstn, aFound, aSucAnt, nClient);

		for (i=0; i<=nClient; i++)
		{
			for (j=0; j<=nClient; j++)
			{
				aTao[i][j]+=aDeltaTaoP[i][j];
				aDeltaTaoP[i][j]=0.0;
			}
		}		
		for (i=0; i<nOptAnt; i++)
		{
			before=0;
			for (k=0; k<=nClient; k++)
			{
				for (j=0; j<=nClient; j++)
				{
					if (aTaboo[optAnt[i]][j]==k)		// �ҵ�һ��·��
					{
						aDeltaTaoP[before][j]=Q3/opt;
						before=j;
						break;
					}
				}
			}
		}
		for (i=0; i<=nClient; i++)
		{
			for (j=0; j<=nClient; j++)
			{
				aTao[i][j]+=aDeltaTaoP[i][j];
				aTao[i][j]=aTao[i][j]*(1.0+(lastopt-opt)/opt);
			}
		}
		lastopt=opt;
	}


	//	������
	nLine=0;
	pResult=fopen("D:\\DataOut.txt", "w");
	if (!pResult)
	{
		return -1;
	}	

	fprintf(pResult, "%d\n", nClient+1);
	for (i=0; i<=nClient; i++)
	{
		fprintf(pResult, "%lf\t%lf\n", aPos[i].x, aPos[i].y);
	}
	fprintf(pResult, "%lf\n", opt);
	cout<<"��̾���ͣ�"<<opt<<endl;
	for (i=0; i<nOptAnt; i++)
	{
//		before=0;
		for (j=1; j<=nClient; j++)
		{
			if (aTaboo[optAnt[i]][j]==0)		// �ҵ�һ��·��
			{
//				fprintf(pResult, "%d-->%d\n", before, j);
				cout<<0<<"-->"<<j;
//				before=j;
				nLine++;
				break;
			}
		}
		for (k=1; k<=nClient; k++)
		{
			for (j=0; j<=nClient; j++)
			{
				if (aTaboo[optAnt[i]][j]==k)		// �ҵ�һ��·��
				{
//					fprintf(pResult, "%d-->%d\n", before, j);
					cout<<"-->"<<j;
//					before=j;
					nLine++;
					break;
				}
			}
		}
		cout<<endl;
	}

	for (i=0; i<NUM_ANT; i++)
	{
		cout<<"����"<<i+1<<"��·��Ϊ��\t"<<0;
		for (k=0; k<=nClient; k++)
		{
			for (j=0; j<=nClient; j++)
			{
				if (aTaboo[i][j]==k)		// �ҵ�һ��·��
				{
//					fprintf(pResult, "%d-->%d\n", before, j);
					cout<<"\t"<<j;
//					before=j;
//					nLine++;
					break;
				}
			}
		}
		cout<<endl;
	}

	fprintf(pResult, "%d\n", nLine);
	for (i=0; i<nOptAnt; i++)
	{
		before=0;
		for (k=0; k<=nClient; k++)
		{
			for (j=0; j<=nClient; j++)
			{
				if (aTaboo[optAnt[i]][j]==k)		// �ҵ�һ��·��
				{
					fprintf(pResult, "%d-->%d\n", before, j);
//					cout<<"-->"<<j;
					before=j;
//					nLine++;
					break;
				}
			}
		}
//		cout<<endl;
	}


	fclose(pResult);
	

	// ����ڴ�
	delete aPos;
	for (i=0; i<=nClient; i++)
	{
		delete aDstn[i];
		delete aDeltaDstn[i];
		delete aTao[i];
		delete aDeltaTaoP[i];
		delete aNTrans[i];
	}
	delete aDstn;
	delete aDeltaDstn;
	delete aTao;
	delete aDeltaTaoP;
	delete aNTrans;
	delete aDmnd;
	delete aCntn;
	for (i=0; i<NUM_ANT; i++)
	{
		delete aTaboo[i];
	}
	delete aTaboo;
	delete aWeight;
	delete aNAnt;
	delete aLastPath;
	delete aProb;
	delete aSuccess;
	delete aFound;
	delete aSucAnt;
	delete optAnt;

	return 0;
}

double Random(double from, double to)
{
	return (double)rand()/(RAND_MAX+1)*(to-from)+from;
}

double FindDistance(Point a, Point b)
{
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}

double FindPower(double in, int power)
{
	int i;
	double result=in;
	if (power<=0)
	{
		return 1;
	}

	for (i=1; i<power; i++)
	{
		result*=in;
	}

	return result;
}

int FindNext(double prob, int nClient, double* aProb/*, bool* aCntn*/)
{
	int i, ret=0;

	for (i=0; i<=nClient; i++)
	{
//		if (aCntn[i])
//		{
		if (prob<aProb[i])
		{
			ret=i;
			break;
		}
//		}
	}

	return ret;
}

// from: �ӵڼ������Ͽ�ʼѰ��·��; nFind: �ҵ��Ŀͻ�����; aTaboo: �����ϵ�·��; aSuccess: �Ƿ�Ϊ����
// Ҫ���·���������ǣ�ֱ���ж���һֻ����; aDeltaTaoP: ����·����ͬʱ����DeltaTaoP; nPath: ����·��
// ��; opt: ���ε����ҵ������·����; aDstn: ����֮���·������; aFound: ���ҵ��Ŀͻ��㼯��; aSucAnt:
// ����Ҫ������ϼ���; nClient: �ͻ�����
void FindSolution(int from, int nFind, const int** aTaboo, bool* aSuccess, double** aDeltaTaoP, int& nPath,
				  double& opt, int* optAnt, int& nOptAnt, int nAnt, double** aDstn, bool* aFound, int* aSucAnt, 
				  const int nClient)
{
	//	double solution;		// �⣨��·�����ȣ�
	int next, i, j, k, before;
	double length;

	while (1)
	{
		while (nFind<nClient)
		{
			next=FindNextAnt(from, nFind, aTaboo, aSuccess, aFound, aSucAnt, nClient, nAnt);
			if (next<0)
			{
				return;
			}

			from=next+1;
			//		FindSolution(next+1, nFind, aTaboo, aSuccess, aDeltaTaoP, nPath, opt, optAnt, 
			//			nOptAnt, nAnt, aDstn, aFound, aSucAnt, nClient);
		}
		//	else
		//	{

		nPath++;
		length=0.0;
		for (i=0; i<nAnt; i++)
		{
			before=0;
			for (k=0; k<=nClient; k++)
			{
				for (j=0; j<=nClient; j++)
				{
					if (aTaboo[aSucAnt[i]][j]==k)		// �ҵ�һ��·��
					{
						length+=aDstn[before][j];
						before=j;
						break;
					}
				}
			}
		}
		for (i=0; i<nAnt; i++)
		{
			before=0;
			for (k=0; k<=nClient; k++)
			{
				for (j=0; j<=nClient; j++)
				{
					if (aTaboo[aSucAnt[i]][j]==k)		// �ҵ�һ��·��
					{
						aDeltaTaoP[before][j]+=Q2/length;
						before=j;
						break;
					}
				}
			}
		}
		if (length<opt)
		{
			opt=length;
			for (i=0; i<nAnt; i++)
			{
				optAnt[i]=aSucAnt[i];
			}
			nOptAnt=nAnt;
		}
		next=FindNextAnt(aSucAnt[nAnt-1]+1, nFind, aTaboo, aSuccess, aFound, aSucAnt, nClient, nAnt);
		if (next<0)
		{
			return;
		}

	}
//	FindSolution(next, nFind, aTaboo, aSuccess, aDeltaTaoP, nPath, opt, optAnt, 
//		nOptAnt, nAnt, aDstn, aFound, aSucAnt, nClient);
	//	}
}

int FindNextAnt(int from, int& nFind, const int** aTaboo, bool* aSuccess, bool* aFound,
				int* aSucAnt, const int nClient, int& nAnt)
{
	int i=from, j, point;
	bool bSuccess=false;

	while (1)
	{
		if (nAnt<0)
		{
			return -1;
		}

		if (nFind>=nClient)
		{
			if (nAnt<=0)
			{
				return -1;
			}
			point=0;
			nAnt--;
			for (j=1; j<=nClient; j++)
			{
				if (aTaboo[aSucAnt[nAnt]][j]>-1)
				{
					aFound[j]=false; 
					point++;
				}
			}
			nFind-=point;

			//		return FindNextAnt(aSucAnt[nAnt]+1, nFind, aTaboo, aSuccess, aFound, aSucAnt, nClient, nAnt);
		}


		while (i<NUM_ANT)
		{
			if (!aSuccess[i])
			{
				i++;
				continue;
			}

			bSuccess=true;
			point=0;
			for (j=1; j<=nClient; j++)
			{
				if (aFound[j]==true && aTaboo[i][j]>-1)
				{
					bSuccess=false;
					break;
				}
				else if (aFound[j]==false && aTaboo[i][j]>-1)
				{
					point++;
				}
			}
			if (!bSuccess)
			{
				i++;
				point=0;
			}
			else
			{
				break;
			}
		}

		if (bSuccess)
		{
			// ����
			for (j=1; j<=nClient; j++)
			{
				if (aTaboo[i][j]>-1)
				{
					aFound[j]=true; 
				}
			}
			aSucAnt[nAnt]=i;
			nAnt++;
			nFind+=point;

			return i;
		}
		else
		{
			// ����
			/*
			if (i<nClient)
			{
			for (j=1; j<=nClient; j++)
			{
			if (aTaboo[i][j]>-1)
			{
			aFound[j]=false; 
			}
			}
			}
			*/
			nAnt--;	
			if (nAnt<0)
			{
				return -1;		// ����
			}
			else
			{		
				point=0;
				for (j=1; j<=nClient; j++)
				{
					if (aTaboo[aSucAnt[nAnt]][j]>-1)
					{
						aFound[j]=false; 
						point++;
					}
				}
				nFind-=point;
				i=aSucAnt[nAnt]+1;
				//			return FindNextAnt(aSucAnt[nAnt]+1, nFind, aTaboo, aSuccess, aFound, aSucAnt, nClient, nAnt);
			}
		}
	}
}