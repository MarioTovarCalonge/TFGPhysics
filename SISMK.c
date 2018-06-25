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
void writeEdgeList(struct WEdge *list[], int ver, int la, int links);
int SFD(double gamma, double a, double b);
void newSFGraph(struct WEdge *list[], int ver, double expo, int la, int *links);

FILE *f;

int link(struct WEdge *list[], int nodei, int nodej);
void update(int N, int L, double *p, double *pav);
void iterate(struct WEdge *list[], int N, int L, double *p, double *pav, double beta, double gamma, double mu);

void writeMatrix(struct WEdge *list[], int ver, int la);

float generador(void);
void iniciar(void);

unsigned int rueda[256];
unsigned char ind_ran,ig1,ig2,ig3;
unsigned int alea;


int main(){
	int t, tmax, vertices, edges, k, l, j, aux;
	double beta, mu, deltabeta, gamma, rhoaux;

	srand(time(NULL));
	iniciar();

    printf("Starting network analisis...\n");
    printf("Enter the Number of Vertices -\n");
    scanf("%d", &vertices);
    l=2;

	double p[vertices*l], pav[vertices*l], rho[l];

	struct WEdge *adlist[l*vertices];

    for (k=0; k<l*vertices; k++){
		adlist[k]=NULL;
    }
	/*Generate or read scale free graph*/

    printf("Reading graph\n");

    double a[2];
    a[0]=2.3;
    a[1]=3.0;

    for(j=0; j<l; j++){

    //readGraph(adlist, vertices, j);
    newSFGraph(adlist, vertices, a[j], j, &edges);
    //writeEdgeList(adlist, vertices, j, edges);
    //writeMatrix(adlist, vertices, j);
    }

    printf("Graph succesfully read\n");


    /*Inicialize network's parameters*/
	tmax=100;
	deltabeta=0.02;
	mu=0.7;
	gamma=0.2;


	f=fopen("resMk.dat", "wt");

	for(beta=0.0; beta<1.0; beta=beta+deltabeta){

        for(k=0; k<vertices*l; k++){
            aux=k/vertices;
            p[k]=0.0;
            pav[k]=0.0;
            rho[aux]=0.0;
        }

        for(k=0; k<(int)(l*vertices*0.05); k++){
            pav[(int)(l*vertices*generador())]=0.0001;
        }

        for(t=0; t<500; t++){
            iterate(adlist, vertices, l, p, pav, beta, gamma, mu);
            update(vertices, l, p, pav);
        }

        rhoaux=0.0;
		for(t=0; t<tmax; t++){
            iterate(adlist, vertices, l, p, pav, beta, gamma, mu);
            update(vertices, l, p, pav);

            for(k=0; k<vertices*l; k++){
                aux=k/vertices;
                rho[aux]=rho[aux]+pav[k]/((double)(tmax));
                rhoaux=rhoaux+pav[k]/((double)(tmax));
                //printf("p[%d]=%lf\n", k, pav[k]);
            }

        }


    if(l==1){
        fprintf(f, "%lf \t %lf \t %lf \t %lf \n", 0.0, 0.0, rhoaux/(l*(double)vertices), beta);
    }
    else{
        fprintf(f, "%lf \t %lf \t %lf \t %lf \n", rho[0]/(vertices), rho[1]/(vertices), rhoaux/(l*(double)vertices), beta);
    }
	printf("|");
	}


fclose(f);
//System is used for calling an extern program for making graphs
//system("plotpantallaMk.gp");
}


struct WEdge *AddWEdge(struct WEdge *currentHead, int newVertex){
    struct WEdge *newHead=(struct WEdge *) malloc(sizeof(struct WEdge));

    newHead->vertex = newVertex;
    newHead->next = currentHead;

    return newHead;
}

//Random graph
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

int SFD(double expo, double a, double b){
    int aux1;
    double aux2;

    do{
    aux1=((int)(b*generador()+a));
    aux2=generador();
    }while(aux2>((1.0)/(pow(aux1, expo))));
return ((int)aux1);
}

//Scale free graph
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

int link(struct WEdge *list[], int nodei, int nodej){
	int r;
	r=0;

	struct WEdge * traverse = list[nodei];

        while (traverse != NULL) {
            if (traverse->vertex==nodej){
				r=1;
				break;
            }
        traverse = traverse->next;
        }
return r;
}

void update(int N, int L, double *p, double *pav){
	int k;

	for(k=0; k<N*L; k++){
		pav[k]=p[k];
		p[k]=0.0;
	}
}

void iterate(struct WEdge *list[], int N, int L, double *p, double *pav, double beta, double gamma, double mu){
	int i, j, l;
	double q;

    for(i=0; i<N*L; i++){
        q=1.0;
        l=i/N;


        struct WEdge * traverse = list[i];

        while (traverse != NULL) {
            q=q*(1.0-beta*pav[traverse->vertex + l*N]);
            traverse = traverse->next;
        }

        if(L!=1){
            if(l==0){
                q=q*(1.0-gamma*pav[i+N]);
            }
            else if(l==L-1){
                q=q*(1.0-gamma*pav[i-N]);
            }
            else{
                q=q*(1.0-gamma*pav[i+N]);
                q=q*(1.0-gamma*pav[i-N]);
            }
        }
        //printf("q=%lf\n", q);

        p[i]=(1.0-q)*(1.0-pav[i])+(1.0-mu)*pav[i]+mu*(1.0-q)*pav[i];
        //p[i]=(1.0-q)+(1.0-mu)*pav[i];
    }

}

void writeMatrix(struct WEdge *list[], int ver, int la){
    int i, j;
    int aux[ver];
    FILE *G;
    char name[100];
    sprintf(name, "A%d.dat", la);
    G=fopen(name, "wt");

    for(i=la*ver;i<(ver+la*ver); i++){
        for(j=0;j<(ver); j++)
            aux[j]=0;

        struct WEdge *traverse=list[i];
        while(traverse!=NULL){
            aux[traverse->vertex]=1;
            traverse=traverse->next;
        }

        for(j=0; j<ver; j++)
            fprintf(G, "%d ", aux[j]);

    fprintf(G, "\n");
    }
fclose(G);
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
