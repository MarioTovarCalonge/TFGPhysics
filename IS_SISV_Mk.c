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
void readGraph(struct WEdge *list[], int N, int links, int la);
int SFD(double gamma, double a, double b);
void newSFGraph(struct WEdge *list[], int ver, double expo, int la, int *links);


FILE *f;
FILE *z;

void update(int N, int L, double *p, double *pav);
void iterate(struct WEdge *list[], int N, int L, double *p, double *pav, double *beta, double gamma, double mu);

float generador(void);
void iniciar(void);

unsigned int rueda[256];
unsigned char ind_ran,ig1,ig2,ig3;
unsigned int alea;


int main(){
	int t, tmax, vertices, edges, k, l, j, aux;
	double mu, gamma;

	srand(time(NULL));
	iniciar();

    printf("Starting network analisis...\n");
    printf("Enter the Number of Vertices -\n");
    scanf("%d", &vertices);
    l=2;


	double p[vertices*l], pav[vertices*l], rho[l];

	struct WEdge *adlist[l*vertices];

	double beta[l];

	beta[1]=0.3;
	beta[0]=beta[1]/6.0;
	gamma=0.25;


    for (k=0; k<l*vertices; k++){
		adlist[k]=NULL;
    }

    /*Generate or read scale free graph*/
    double a[2];
    a[0]=2.7;
    a[1]=2.7;
    printf("Reading graph\n");

    for(j=0; j<l; j++){
    //newGraph(Wadlist, vertices, edges, "undirected", j);
    newSFGraph(adlist, vertices, a[j], j, &edges);
    //writeEdgeList(Wadlist, vertices, j);
    }
    printf("Graph succesfully read\n");


    /*Inicialize network's parameters*/
	tmax=50;
	mu=0.25;

	f=fopen("capa2Mk.dat", "wt");
	z=fopen("capa1Mk.dat", "wt");

        for(k=0; k<vertices*l; k++){
            aux=k/vertices;
            p[k]=0.0;
            pav[k]=0.0;
            rho[aux]=0.0;
        }

        for(k=0; k<(int)(vertices*0.1); k++){
            aux=(int)(vertices*generador()+vertices);
            pav[aux]=1.0;
            pav[aux-vertices]=1.0;
        }


		for(t=0; t<tmax; t++){

            rho[0]=0.0;
            rho[1]=0.0;

            iterate(adlist, vertices, l, p, pav, beta, gamma, mu);
            update(vertices, l, p, pav);

            for(k=0; k<vertices*l; k++){
                aux=k/vertices;
                rho[aux]=rho[aux]+pav[k];
                //printf("p[%d]=%lf\n", k, pav[k]);
            }

        fprintf(f, "%lf \t %lf \t %d\n", rho[1]/((double)vertices), 1.0-rho[1]/((double)vertices), t);
        fprintf(z, "%lf \t %lf \t %d \n", rho[0]/((double)vertices), 1.0-rho[0]/((double)vertices), t);
        printf("|");
        }



fclose(f);
fclose(z);
system("plotMK.gp");
}


struct WEdge *AddWEdge(struct WEdge *currentHead, int newVertex){
    struct WEdge *newHead=(struct WEdge *) malloc(sizeof(struct WEdge));

    newHead->vertex = newVertex;
    newHead->next = currentHead;

    return newHead;
}

void newGraph(struct WEdge *list[], int ver, int enlaces, char *dir, int la){
	int aux1, aux2, i, r;

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

void readGraph(struct WEdge *list[], int N, int links, int la){
    int i, v1, v2;
    FILE *g;
    char aux[100];
    sprintf(aux, "EdgeList%d.dat", la+1);
    g=fopen(aux, "rt");
    for(i=0; i<links; i++){
        fscanf(g, "%d%d", &v1, &v2);
        list[v1+la*N] = AddWEdge(list[v1+la*N], v2);
    }

fclose(g);
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

void update(int N, int L, double *p, double *pav){
	int k;

	for(k=0; k<N*L; k++){
		pav[k]=p[k];
		p[k]=0.0;
	}
}

void iterate(struct WEdge *list[], int N, int L, double *p, double *pav, double *beta, double gamma, double mu){
	int i, l;
	double q;

    for(i=0; i<N*L; i++){
        q=1.0;
        l=i/N;

        struct WEdge * traverse = list[i];

        if(i<N){
            while (traverse != NULL) {
                q=q*(1.0-beta[0]*pav[traverse->vertex]);
                traverse = traverse->next;
            }
            q=q*(1.0-1*pav[i+N]);

            p[i]=(1.0-q)*(1.0-pav[i])+pav[i];
        }
        else{
            while (traverse != NULL) {
                q=q*(1.0-beta[1]*pav[traverse->vertex + N]);
                traverse = traverse->next;
            }
            //q=q*gamma*pav[i-N];
            //if(q<0.0)
                //q=0.0;

            //printf("q=%lf\n", q);

            p[i]=(1.0-q)*(1.0-pav[i])+(1.0-mu)*pav[i]-gamma*pav[i-N];//+mu*(1.0-q)*pav[i]
            if(p[i]<0.0)
                p[i]=0.0;
            //printf("p=%lf\n", p[i]);
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

float generador(void){
	float r;

	ig1=ind_ran-24;
	ig2=ind_ran-55;
	ig3=ind_ran-61;
	rueda[ind_ran]=rueda[ig1]+rueda[ig2];
	alea=(rueda[ind_ran]^rueda[ig3]);
	ind_ran++;

	r=alea*NormRANu;
	return r;
}
