#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SQR(x) ((x)*(x)) //square

struct Nombre_config
{
  int n_in;        // NOMBRE_DE_ATTRIBUTS
  int n_l_out;     // NOMBRE_DE_LIGNES_DES_NEURONES
  int n_c_out;     // NUMBER_DE_COLONNES_DES_NEURONES
  int n_out;       // NUMBER_TOTAL_DE_NEURONES (n_l_out * n_c_out)
  int nb_it;       // NOMBRE_ITERATION
  double minAlpha; // ALPHA_PHASE (valeur de départ)
  int Etti;       // Ettiquettage
  int pEtti;
}Nb_config;


struct node  //neuron (node)
{
  double act;   // LA DISTANCE EUCLIDIENNE
  char *etiq;   // ETTIQUETTAGE
  double *w;   //  POIDS DU VECTEUR
};

typedef struct node t_node;

struct bmu {    // bmu Coordonnees
double act; // LA DISTANCE EUCLIDIENNE
int r;
int c;
};

typedef struct bmu t_bmu;  // cherche Le BMU
t_bmu *Bmu = NULL;
int Bmu_size=1;



struct carte   //
{
  int rayon_voisinage;  // rayon_voisinage
  t_node **map;      // matrice
  double *d_vect;  //  data vector
  double alpha;     // coeficient d'apprentissage
  char *etiq;       // ETTIQUETTAGE
} Carte;


struct vecteur
{
        double *mat;
        char *name;
        double norm;
};

struct vecteur * matrice_vecteur;




// STRUCT LES VECTEURS DE LA MATRICE
void struct_vecteur(int n)
{
    matrice_vecteur=malloc(n*sizeof(struct vecteur));
    int i;
    for(i=0;i<n;i++)
    {
        matrice_vecteur[i].mat=malloc(Nb_config.n_in*sizeof(double));
        matrice_vecteur[i].name=malloc(20*sizeof(char));
    }
}



double *moyenne,*min,*max;

//UNE FONCTION POUR PRODUIRE LE VECTEUR MOYEN A PARTIR D"UNE MATRICE
void moyen_De_Vecteur(int n)
{
    moyenne=malloc(Nb_config.n_in*sizeof(double));
    memset(moyenne,0,Nb_config.n_in*sizeof(double));

    int i,j;

    for(i=0;i<Nb_config.n_in;i++)
    {
        for(j=0;j<n;j++)
            moyenne[i]+=matrice_vecteur[j].mat[i];
        moyenne[i]/=n;
    }
}



//UNE FONCTION POUR PRODUIRE valeur Min
void min_vec(double k)
{
    min=malloc(Nb_config.n_in*sizeof(double));
    int i;
    for(i=0;i<Nb_config.n_in;i++)
        min[i]=moyenne[i]-k;
}



//UNE FONCTION POUR PRODUIRE valeur Max
void max_vec(double k)
{
    max=malloc(Nb_config.n_in*sizeof(double));
    int i;
    for(i=0;i<Nb_config.n_in;i++)
        max[i]=moyenne[i]+k;
}


//UNE FONCTION POUR NORMALISER UNE MATRICE
void normelise_matrice_vecteur(int i,int size)
{
    double sum=0.;
    int j;
    for(j=0;j<size;j++)
        sum+=SQR(matrice_vecteur[i].mat[j]);
    matrice_vecteur[i].norm=sqrt(sum);
}

//UNE FONCTION POUR DENORMALISER UNE MATRICE
void denormalise_matrice_vecteur(int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<Nb_config.n_in;j++)
            matrice_vecteur[i].mat[j]/=matrice_vecteur[i].norm;
    }
}

double* init_rand_w()
{
    int i;
    double k=(double)rand()/RAND_MAX;
    double *tmp_w=malloc(Nb_config.n_in*sizeof(double));

    for(i=0;i<Nb_config.n_in;i++)
        {
            tmp_w[i]=k*(max[i]-min[i])+min[i];
        }

    double norm=0.;

    for(i=0;i<Nb_config.n_in;i++)
    {
        norm+=SQR(tmp_w[i]);
    }

    for(i=0;i<Nb_config.n_in;i++)
    {
            tmp_w[i]/=norm;
    }
    return tmp_w;
}


// MÉLANGER MON ESPACE DE DONNÉES  :

int * index_matrice;
void melangeToi(int n)
{
    index_matrice=malloc(sizeof(int)*n);
    int i;
    for(i=0;i<n;i++)
        index_matrice[i]=i;
}



void melange_matrice(int n)
{
    int i,r_and,k;
    srand(time(NULL));
    for(i=0;i<n;i++)
        {
            r_and=rand()%n;
            k=index_matrice[i];
            index_matrice[i]=index_matrice[r_and];
            index_matrice[r_and]=k;
        }
}



 //UNE FONCTION QUI CALCULE LA DISTANCE EUCLIDIENNE ENTRE DEUX VECTEURS
double euc_distance(double *a1, double *a2, int n)
{
double sum=0.;
int i;
for(i=0;i<n;i++)
{
sum+=(SQR(a1[i] - a2[i]));
}
return sqrt(sum);
}



void calc_alpha(int it_n, int tot_it)
{
Carte.alpha = Nb_config.minAlpha * (1 - ((double)it_n/(double)tot_it));
}



void update(t_bmu* b_mu)
{
    int nr=Carte.rayon_voisinage;
    int i,j,x1,x2,y1,y2;//top and bottom

    for(;nr>=0;nr--)
    {
        if(b_mu->r-nr<0)
            x1=0;
        else
            x1=b_mu->r-nr;
        if(b_mu->c-nr<0)
            y1=0;
        else
            y1=b_mu->c-nr;

        if(b_mu->r+nr>Nb_config.n_l_out-1)
            x2=Nb_config.n_l_out-1;
        else
            x2=b_mu->r+nr;
        if(b_mu->c+nr>Nb_config.n_c_out-1)
            y2=Nb_config.n_c_out-1;
        else
            y2=b_mu->c+nr;

        for(i=x1;i<=x2;i++)
            for(j=y1;j<=y2;j++)
            {

                int k;

                for(k=0;k<Nb_config.n_in;k++)
                    {

                        Carte.map[i][j].w[k]+=Carte.alpha*(Carte.d_vect[k]-Carte.map[i][j].w[k]);
                    }
            }
    }
}


// INITIALISATION LES NOMBRES LA DATA
void init_Nombre_config()
{
 Nb_config.n_l_out=6;
 Nb_config.n_c_out=10;
 Nb_config.n_out=Nb_config.n_l_out*Nb_config.n_c_out;
 Nb_config.n_in=4;
 Nb_config.nb_it=30000;
 Nb_config.minAlpha=0.7;
 Nb_config.pEtti=Nb_config.nb_it/5;
 Nb_config.Etti=2;
}



// UNE FONCTION QUI LIT LE CONTENU D'UN FICHIER
void lecteur_data()
{
    FILE * in;

    char *str=malloc(sizeof(char)*500);
    in=fopen("iris.data","r");

    int i,j;
for(i=0;i<150;i++)
{
        fscanf(in,"%s",str);
        char *tok=strtok(str,",");

        for(j=0;j<Nb_config.n_in;j++)
            {
                matrice_vecteur[i].mat[j]=atof(tok);
                tok=strtok(NULL,",");
            }

        if (strcmp(tok, "Iris-setosa") == 0)
strcpy(matrice_vecteur[i].name,"A");
        else if(strcmp(tok,"Iris-versicolor")==0)
            strcpy(matrice_vecteur[i].name,"B");
        else
            strcpy(matrice_vecteur[i].name,"C");

        normelise_matrice_vecteur(i,Nb_config.n_in);
}

fclose(in);
    free(str);
}


// CREATION LES NEURON DE LA MATRICE
void create_neuron_matrice()
{
    int i,j;
    Carte.map=malloc(Nb_config.n_l_out*sizeof(t_node *));
for(i=0;i<Nb_config.n_l_out;i++)
{
Carte.map[i]=malloc(Nb_config.n_c_out*sizeof(t_node));
}
for(i=0;i<Nb_config.n_l_out;i++)
{
for (j=0;j<Nb_config.n_c_out;j++)
{

            Carte.map[i][j].w=(double*)malloc(sizeof(double)*Nb_config.n_in);
Carte.map[i][j].w=init_rand_w();
Carte.map[i][j].etiq=malloc(20*sizeof(char));
strcpy(Carte.map[i][j].etiq, ".");
}
}
}


// FONCTION POUR AFFICHE LE CONTENU DE LA MATRICE

void AfficheContenuDeMatrice()
{
    int i,j;
    for(i=0;i<Nb_config.n_l_out;i++)
    {
        for(j=0;j<Nb_config.n_c_out;j++)
            {
                printf("%s ",Carte.map[i][j].etiq);
            }
        printf("\n");
    }
}




void Ettiquettage()
{
    int i,j,p,u,it;
    double min_d,dist;



    Bmu=malloc(sizeof(t_bmu));

    for(p=0;p<Nb_config.Etti;p++)
    {
        int cur_n_it;
        if(!p)
        {
            cur_n_it=Nb_config.pEtti;
        }
        else
        {
            cur_n_it=Nb_config.nb_it-Nb_config.pEtti;
            Nb_config.minAlpha=0.07;
            Carte.rayon_voisinage=1;
        }

        for(it=0;it<cur_n_it;it++)
        {
            calc_alpha(it,cur_n_it);

            if(it%(Nb_config.pEtti/2)==0&&it!=0&&p==0)
            {
                Carte.rayon_voisinage-=1;
            }


            melange_matrice(150);

            for(u=0;u<150;u++)
            {

                Carte.d_vect=matrice_vecteur[index_matrice[u]].mat;
                min_d=1000.;
                for(i=0;i<Nb_config.n_l_out;i++)
                {
                    for(j=0;j<Nb_config.n_c_out;j++)
                    {
                        dist=euc_distance(Carte.d_vect,Carte.map[i][j].w,Nb_config.n_in);
                        Carte.map[i][j].act=dist;
                        if(dist<min_d)
                        {
                            min_d=dist;
                            if(Bmu_size>1)
                            {
                                Bmu_size=1;
                                Bmu=realloc(Bmu,Bmu_size*sizeof(t_bmu));
                            }
                            Bmu[0].act=dist;
                            Bmu[0].r=i;
                            Bmu[0].c=j;
                        }
                        else if(dist==min_d)
                        {

                            Bmu_size++;
                            Bmu=realloc(Bmu,Bmu_size*sizeof(t_bmu));
                            Bmu[Bmu_size-1].act=dist;
                            Bmu[Bmu_size-1].r=i;
                            Bmu[Bmu_size-1].c=j;

                        }
                    }
                }

                if(Bmu_size>1)
                {
                    int t=rand()%(Bmu_size);
                    Bmu[0]=Bmu[t];
                }

                strcpy(Carte.map[Bmu[0].r][Bmu[0].c].etiq, matrice_vecteur[index_matrice[u]].name);
                update(Bmu);
            }
        }
    }
}




int main()
{
    init_Nombre_config();

    struct_vecteur(150);

lecteur_data();
    denormalise_matrice_vecteur(150);

    moyen_De_Vecteur(150);
    min_vec(0.004);
    max_vec(0.006);

    melangeToi(150);

    // créer une matrice de neurones avec un neurone aléatoire

create_neuron_matrice();
printf("Affiche Contenu de la matrice d'ettiquettage par des poins:\n");
AfficheContenuDeMatrice();
    printf("\n");

    Carte.rayon_voisinage=6;
    Carte.alpha=0;

    Ettiquettage();
    printf("Affiche Le Contenu final de la matrice :\n");
    AfficheContenuDeMatrice();
    free(moyenne);
    free(min);
    free(max);

return 0;
}
