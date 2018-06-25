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

struct promedio{
    double sano;
    double vacunado;
    double infectado;
};


struct WEdge *AddWEdge(struct WEdge *currentHead, int newVertex);


void newGraph(struct WEdge *list[], int ver, int enlaces, char *dir, int la);
void readGraph(struct WEdge *list[], int links);
void writeEdgeList(struct WEdge *list[], int ver, int la);
int SFD(double gamma, double a, double b);
void newSFGraph(struct WEdge *list[], int ver, double expo, int la, int *links);


FILE *f;
FILE *z;

void update(int N, int L, int *p, int *pav);
void iterateMC(struct WEdge *list[], int N, int L, int *p, int *pav, double *beta, double *gamma, double mu, double omega);


double generador(void);
void iniciar(void);

unsigned int rueda[256];
unsigned char ind_ran,ig1,ig2,ig3;
unsigned int alea;


int main(){
	int vertices, edges, k, l, j, rep, repmax, MC, MCmax, aux;
	double mu, omega;
	int est[3][2];

	srand(time(NULL));
	iniciar();

    printf("Starting network analisis...\n");
    printf("Enter the Number of Vertices -\n");
    scanf("%d", &vertices);
	l=2;

	int p[vertices*l], pav[vertices*l];
	struct WEdge *Wadlist[l*vertices];
	double gamma[l], beta[l];

	beta[1]=0.6;
	beta[0]=0.4;
	//beta[0]=0;

	gamma[1]=1.0;
	gamma[0]=0.15;

    for (k=0; k<vertices*l; k++){
		Wadlist[k]=NULL;
    }

	/*Generate or read scale free graph*/
    double a[2];
    a[0]=2.7;
    a[1]=2.7;
    printf("Reading graph\n");
    for(j=0; j<l; j++){
    //newGraph(Wadlist, vertices, edges, "undirected", j);
    newSFGraph(Wadlist, vertices, a[j], j, &edges);
    //writeEdgeList(Wadlist, vertices, j);
    }
    printf("Graph succesfully read\n");

	/*Inicialize network's parameters*/
	mu=0.3;
	omega=0.15;

	repmax=20;
    MCmax=50;

    struct promedio med[MCmax*l];

    for(j=0; j<MCmax*l; j++){
        med[j].sano=0.0;
        med[j].vacunado=0.0;
        med[j].infectado=0.0;
    }


    f=fopen("capa2.dat", "wt");
    z=fopen("capa1.dat", "wt");

		for(rep=0; rep<repmax; rep++){

            /*Inicialize probability arrays*/
            for(k=0; k<l*vertices; k++){
                p[k]=0;
                pav[k]=0;
            }

			//Random start of infected nodes, five percent of total network
			for(k=0; k<(int)(vertices*0.1); k++){
                aux=(int)(vertices*generador()+vertices);
				pav[aux]=1;
                pav[aux-vertices]=1;

			}


			for (MC=0; MC<MCmax; MC++){
                for(k=0; k<3; k++){
                    est[k][0]=0;
                    est[k][1]=0;
                }

				iterateMC(Wadlist, vertices, l, p, pav, beta, gamma, mu, omega);
				update(vertices, l, p, pav);

				for(k=0; k<l*vertices; k++){
				    aux=k/vertices;
                    if(pav[k]==1){
                        est[1][aux]++;
                    }
                    if(pav[k]==3){
                        est[2][aux]++;
                    }
                    if(pav[k]==0){
                        est[0][aux]++;
                    }
				}

                med[MC].sano+=est[0][0]/((double)vertices);
                med[MC].vacunado+=est[2][0]/((double)vertices);
                med[MC].infectado+=est[1][0]/((double)vertices);
                med[MC+MCmax].sano+=est[0][1]/((double)vertices);
                med[MC+MCmax].vacunado+=est[2][1]/((double)vertices);
                med[MC+MCmax].infectado+=est[1][1]/((double)vertices);

			}
        printf("|");
		}

		for (MC=0; MC<MCmax; MC++){
            fprintf(f, "%lf \t %lf \t %lf \t %d \t %lf\n", med[MC+MCmax].sano/((double)repmax), med[MC+MCmax].infectado/((double)repmax), med[MC+MCmax].vacunado/((double)repmax), MC, med[MC+MCmax].sano/repmax+med[MC+MCmax].infectado/repmax+med[MC+MCmax].vacunado/repmax);
            fprintf(z, "%lf \t %lf \t %lf \t %d \t %lf\n", med[MC].sano/((double)repmax), med[MC].infectado/((double)repmax), med[MC].vacunado/((double)repmax), MC, med[MC].sano/repmax+med[MC].infectado/repmax+med[MC].vacunado/repmax);
		}

fclose(f);
fclose(z);
system("plotMC.gp");
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

void readGraph(struct WEdge *list[], int links){
    int i, v1, v2;
    FILE *g;
    g=fopen("Input.dat", "rt");
    for(i=0; i<links; i++){
        fscanf(g, "%d%d%d", &v1, &v2);
        list[v1] = AddWEdge(list[v1], v2);
    }

fclose(g);
}

void writeEdgeList(struct WEdge *list[], int ver, int la){
    int i;
    struct WEdge *go;
    FILE *H;
    char aux[100];
    sprintf(aux, "EdgeList%d.dat", la+1);
    H=fopen(aux, "wt");

    for(i=0+la*ver; i<ver+la*ver; i++){
        go=list[i];
        while(go!=NULL){
            fprintf(H, "%d %d\n", i-la*ver, go->vertex);
            go=go->next;
        }
    }
fclose(H);
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

void update(int N, int L, int *p, int *pav){
	int i;

	for(i=0; i<N*L; i++){
		pav[i]=p[i];
		p[i]=0;
	}
}

void iterateMC(struct WEdge *list[], int N, int L, int *p, int *pav, double *beta, double *gamma, double mu, double omega){
	int k, aux;

		for(k=0; k<N*L; k++){
            aux=k/N;
			if(pav[k]==1){

				struct WEdge * traverse = list[k];

				while (traverse != NULL){

                    if (generador()<beta[aux]){
                        if(pav[traverse->vertex + aux*N]==0){
                            p[traverse->vertex + aux*N]=1;
                        }
					}

                traverse = traverse->next;
				}

                if(k<N){
                    if(generador()<omega){
                        p[k]=3;
                    }
					else{
						p[k]=1;
					}

                    if(generador()<gamma[0]){
                        pav[k+N]=3;
                    }
                }
                else{
					if(pav[k-N]!=3)
                        p[k-N]=1;
					
                    if(generador()<mu){
                        p[k]=0;
                    }
					else{
						p[k]=1;
					}
                }

			}
			else{
                if(pav[k]==3)
                    p[k]=3;
			}
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
