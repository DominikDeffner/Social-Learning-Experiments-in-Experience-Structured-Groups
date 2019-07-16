// Social Learning with frequency and experience bias

data{
int N;
int N_id;
int Session_ID[N];
int id[N];
int trial[N];
int group_id[N];
int Choice[N];
int Correct[N];
real Payoff[N];
int Experience[N];
int farm[N];
int new_farm[N];
int n1[N];
int n2[N];
int n3[N];
int n4[N];

matrix[N,4] nmat;

// age vars

matrix[N,3] age_models;
matrix[N,3] choice_models;
}

parameters{
real phi;
real<lower=0> L;
real<lower=0,upper=1> sigma; // weight of social learning
real<lower=0> f; // strength of conformity
real<lower=0,upper=1> kappa; // weight of age bias
real beta; // strength of age bias

}

model{

matrix[N_id,4] A; // attraction matrix

phi ~ beta(2,2);
L ~ exponential(1);
sigma ~ beta(2,2);
f ~ normal(1,0.5);
kappa ~ beta(2,2);
beta ~ normal(0,0.5);

// initialize attraction scores

for ( i in 1:N_id ) A[i,1:4] = rep_vector(0,4)';

// loop over choices

for ( i in 1:N ) {

vector[4] pay;
vector[4] pA;
vector[4] pC;
vector[4] pS;
vector[4] p;

if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

// first, what is log-prob of observed choice

pA = softmax( L*A[id[i],1:4]' );

// second compute conformity thing

if ( sum(nmat[i,:])==0 ) {
p = pA;

} else {

// conformity

for ( j in 1:4 ) pC[j] = nmat[i,j]^f;

pC = pC / sum(pC);

//age bias

for ( an_option in 1:4 ){

pS[an_option] = 0;

for ( a_model in 1:3 ) {

if ( choice_models[i,a_model]==an_option )
pS[an_option] = pS[an_option] + exp(beta*age_models[i,a_model]);

}
}

pS = pS / sum(pS);

// combine everything

p = (1-sigma)*pA + sigma*( (1-kappa)*pC + kappa*pS );

}

Choice[i] ~ categorical( p );

// second, update attractions conditional on observed choice

pay[1:4] = rep_vector(0,4);
pay[ Choice[i] ] = Payoff[i];
A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';


}
//i
}
