#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define NormRANu (2.3283063671E-10F)


struct WEdge {
    int vertex;
    struct WEdge *next;
};


struct WEdge *AddWEdge(struct WEdge *currentHead, int newVertex);


void newGraph(struct WEdge *list[], int ver, int enlaces, char *dir, int la);
void readGraph(struct WEdge *list[], int N, int la);
int SFD(double gamma, double a, double b);
void newSFGraph(struct WEdge *list[], int ver, double expo, int la, int *links);
void writeEdgeList(struct WEdge *list[], int ver, int la, int links);


FILE *f;
FILE *z;

void update(int N, int L, int *p, int *pav);
void iterateMC(struct WEdge *list[], int N, int L, int *p, int *pav, double beta, double gamma, double mu);


double generador(void);
void iniciar(void);

unsigned int rueda[256];
unsigned char ind_ran,ig1,ig2,ig3;
unsigned int alea;


int main(){
	int vertices, edges, k, l, j, rep, repmax, term, MC, MCmax, media, aux1;
	double beta, mu, deltabeta, gamma, var, tol;
	double rhoR, rhoMC, rho1, rho2, rhoR1, rhoR2;

	srand(time(NULL));
	iniciar();

    printf("Starting network analisis...\n");
    printf("Enter the Number of Vertices -\n");
    scanf("%d", &vertices);
	l=2;

	int p[vertices*l], pav[vertices*l];
	int prob[vertices*l];
	struct WEdge *Wadlist[l*vertices];

    // Hay que inicializar el puntero de la lista
    for (k=0; k<vertices*l; k++){
		Wadlist[k]=NULL;
    }

	/*Generate or read scale free graph*/

    printf("Reading graph\n");
    double a[2];
    a[0]=2.3;
    a[1]=3.0;

    for(j=0; j<l; j++){
    //newGraph(Wadlist, vertices, edges, "undirected", j);
    newSFGraph(Wadlist, vertices, a[j], j, &edges);
    //writeEdgeList(Wadlist, vertices, j, edges);
    //readGraph(Wadlist, vertices, j);
    }
    printf("Graph succesfully read\n");

	/*Inicialize network's parameters*/
	deltabeta=0.02;
	mu=0.7;
	gamma=0.2;
	repmax=50;
	tol=0.00001;

	//El indice del nodo nos da el link interlayer, intralayer se propaga como en singlelayer.

	f=fopen("resMC.dat", "wt");
	z=fopen("datarep.dat", "wt");

    MCmax=50;
    double mediaRep[repmax];
    for(k=0; k<repmax; k++){
        mediaRep[k]=0.0;
    }
	//We vary beta for studying all possibilities
	for(beta=0.0; beta<1.0; beta=beta+deltabeta){
		//We repeat for statistic purposes
		rhoR=0.0;
		rhoR1=0.0;
		rhoR2=0.0;
		media=0.0;


		for(rep=0; rep<repmax; rep++){
            //g=fopen("dataMC.dat", "wt");
            /*Inicialize probability arrays*/
            for(k=0; k<l*vertices; k++){
                p[k]=0;
                pav[k]=0;
                prob[k]=0;
            }

			//Random start of infected nodes, five percent of total network
			for(k=0; k<(int)(l*vertices*0.05); k++){
				pav[(int)(l*vertices*generador())]=1;
			}

			for(term=0; term<200; term++){
                iterateMC(Wadlist, vertices, l, p, pav, beta, gamma, mu);
				update(vertices, l, p, pav);
			}


			for (MC=0; MC<MCmax; MC++){
				iterateMC(Wadlist, vertices, l, p, pav, beta, gamma, mu);
				update(vertices, l, p, pav);

                for(k=0; k<l*vertices; k++){
                    if(pav[k]!=0)
                        prob[k]++;
                }
			}
            rhoMC=0.0;
            rho1=rho2=0.0;
			for(k=0; k<l*vertices; k++){
                rhoMC=rhoMC+prob[k]/((double)MCmax);
                if(k<vertices){
                    rho1=rho1+prob[k]/((double)MCmax);
                }
                else{
                    rho2=rho2+prob[k]/((double)MCmax);
                }
            }

		fprintf(z, "%lf \t %lf \n", rhoMC/((double)(l*vertices)), beta);

		if((rhoMC/((double)(l*vertices)))>tol){
            rhoR=rhoR+rhoMC/((double)(l*vertices));
            rhoR1=rhoR1+rho1/((double)(vertices));
            rhoR2=rhoR2+rho2/((double)(vertices));
            mediaRep[rep]=rhoMC/((double)(l*vertices));
            media++;
		}
		else{
            if(media!=0.0)
                rep=rep-1;
		}

    }
    var=0.0;
    aux1=0;
    for(j=0; j<repmax; j++){
        if(mediaRep[j]>0.0){
            var=var+(mediaRep[j]-rhoR/((double)(media)))*(mediaRep[j]-rhoR/((double)(media)));
            aux1++;
        }
    }
    if(aux1!=0)
        var=var/((double)aux1);

    if(media==0)
        media++;

	fprintf(f, "%lf \t %lf \t %lf \t %lf \t %.10lf\n", rhoR/((double)(media)), rhoR1/((double)(media)), rhoR2/((double)(media)),  beta, var);
	printf("|");
	}


fclose(f);
fclose(z);
//fclose(g);
system("plotpantalla.gp");
}


struct WEdge *AddWEdge(struct WEdge *currentHead, int newVertex){
    struct WEdge *newHead=(struct WEdge *) malloc(sizeof(struct WEdge));

    newHead->vertex = newVertex;
    newHead->next = currentHead;

    return newHead;
}

void newGraph(struct WEdge *list[], int ver, int enlaces, char *dir, int la){
	int aux1, aux2, i, r;

		iniciar();

        for(i=0; i<enlaces; i++){
			//Generar vertices entre los que crear enlaces
            aux1=(int)(ver*generador());
            do{
                aux2=(int)(ver*generador());
            }while(aux1==aux2); //Evitamos self-edges

            r=0;
            struct WEdge * go = list[aux1+la*ver];
            while(go!=NULL){
                if(go->vertex==aux2)
                    r=1;
            go = go->next;
            }

            if(r==0){
                list[aux1+la*ver] = AddWEdge(list[aux1+la*ver], aux2);

                if(strcmp(dir, "undirected")==0)
                    list[aux2+la*ver] = AddWEdge(list[aux2+la*ver], aux1);
            }
            else{
                i=i-1;
            }
        }

}

int SFD(double expo, double a, double b){
    int aux1;
    double aux2;

    do{
    aux1=((int)(b*generador()+a));
    aux2=generador();
    }while(aux2>((1.0)/(pow(aux1, expo))));
return ((int)aux1);
}

void newSFGraph(struct WEdge *list[], int ver, double expo, int la, int *links){
    int aux[ver];
    int i, k, aux1, aux2, r, enlaces;
    enlaces=0;

    for(i=0; i<ver; i++){
        aux[i]=SFD(expo, 1, ver);
        if(aux[i]/2<aux[i]/2.0){
            if((aux[i]+1)<ver){
                aux[i]=aux[i]+1;
            }
            else{
                aux[i]=aux[i]-1;
            }

        }
        enlaces=enlaces+aux[i];
    }
        *links=enlaces;

        for(i=0; i<ver; i++){
			for(k=0; k<aux[i]; k++){

            aux1=i;
            do{
                aux2=(int)(ver*generador());
            }while(aux1==aux2); //Evitamos self-edges

            r=0;
            struct WEdge * go = list[aux1+la*ver];
            while(go!=NULL){
                if(go->vertex==aux2)
                    r=1;
            go = go->next;
            }

            if(r==0){
                list[aux1+la*ver] = AddWEdge(list[aux1+la*ver], aux2);
                list[aux2+la*ver] = AddWEdge(list[aux2+la*ver], aux1);
            }
            else{
                k=k-1;
            }
            }
        }
}

void readGraph(struct WEdge *list[], int N, int la){
    int i, v1, v2, enlaces;
    FILE *g;
    char aux[100];
    sprintf(aux, "EdgeList%d.dat", la+1);
    g=fopen(aux, "rt");
    fscanf(g, "%d%d", &v1, &v2);
    enlaces=v2;
    for(i=0; i<enlaces*2; i++){
        fscanf(g, "%d%d", &v1, &v2);
        list[v1+la*N] = AddWEdge(list[v1+la*N], v2);
    }

fclose(g);
}

void writeEdgeList(struct WEdge *list[], int ver, int la, int links){
    int i;
    struct WEdge *go;
    FILE *H;
    char aux[100];
    sprintf(aux, "EdgeList%d.dat", la+1);
    H=fopen(aux, "wt");
    fprintf(H, "%d %d\n", ver, links);

    for(i=0+la*ver; i<ver+la*ver; i++){
        go=list[i];
        while(go!=NULL){
            fprintf(H, "%d %d\n", i-la*ver, go->vertex);
            go=go->next;
        }
    }
fclose(H);
}

void update(int N, int L, int *p, int *pav){
	int i;

	for(i=0; i<N*L; i++){
		pav[i]=p[i];
		p[i]=0;
	}
}

void iterateMC(struct WEdge *list[], int N, int L, int *p, int *pav, double beta, double gamma, double mu){
	int k, control;

		for(k=0; k<N*L; k++){

            control=k/N;

		//MC iteration for node k
			if(pav[k]!=0){
				//We try to infect the neibourghs of node k
				struct WEdge * traverse = list[k];//list depends of layer under study

				while (traverse != NULL) {
                    if (generador()<beta){
                        if(pav[traverse->vertex + control*N]!=0){
						p[traverse->vertex + control*N]=2;
                        }
                        else{
                        p[traverse->vertex + control*N]=1;
                        }
					}
                    traverse = traverse->next;
				}

				//We try to infect node's counterpart in other layers
                if(L!=1){
                if(control==0){
                    if (generador()<gamma){
                        if(p[k + N]==0){
                            p[k + N]=1;
                        }
                        if(p[k + N]==1){
                            p[k + N]=2;
                        }
                    }
                }
                else if(control==L-1){
                    if (generador()<gamma){
                        if(p[k - N]==0){
                            p[k - N]=1;
                        }
                        if(p[k - N]==1){
                            p[k - N]=2;
                        }
                    }
				}
				else{
                    if (generador()<gamma)
                        p[k + N]=1;
                    if (generador()<gamma)
                        p[k - N]=1;
				}
                }

				if (generador()<mu){
                    if(p[k]==2){
                        p[k]=1;
                    }
                    else{
                        p[k]=0;
                    }
				}
				else{
                    p[k]=pav[k];
				}

			}
        //printf(" iteracion %d bien ", k);
		}

}

void iniciar(void){
    int i;

    for(i=0; i<256; i++){
		rueda[i]=((rand()<<16)+rand());
	}
    ind_ran=0;
    ig1=0;
    ig2=0;
    ig3=0;
}

double generador(void){
	double r;

	ig1=ind_ran-24;
	ig2=ind_ran-55;
	ig3=ind_ran-61;
	rueda[ind_ran]=rueda[ig1]+rueda[ig2];
	alea=(rueda[ind_ran]^rueda[ig3]);
	ind_ran++;

	r=alea*NormRANu;
	return r;
}
