/**
	@file
	mda_stereo~: a simple audio object for Max
	original by: jeremy bernstein, jeremy@bootsquad.com
	@ingroup examples
*/

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "ext_buffer.h"
#include "z_dsp.h"			// required for MSP objects




// struct to represent the object's state
typedef struct _mda_stereo {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
    t_double fli, fld, fri, frd, fdel, phi, dphi, mod;
    t_double fParam1, fParam2, fParam3, fParam4, fParam5; // width, delay, balance, mod, rate
    t_symbol *mda_stereo_mode;
    
    t_int32 size, bufpos;
    t_double *buffer;
} t_mda_stereo;


// method prototypes
void *mda_stereo_new(t_symbol *s, long argc, t_atom *argv);
void mda_stereo_free(t_mda_stereo *x);
void mda_stereo_assist(t_mda_stereo *x, void *b, long m, long a, char *s);
void mda_stereo_width(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv);
void mda_stereo_delay(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv);
void mda_stereo_balance(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv);
void mda_stereo_mod(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv);
void mda_stereo_rate(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv);
void mda_stereo_mode(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv);
void mda_stereo_reset(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv);
void mda_stereo_recalculate(t_mda_stereo *x);
void mda_stereo_dsp64(t_mda_stereo *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void mda_stereo_perform64(t_mda_stereo *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *mda_stereo_class = NULL;


//***********************************************************************************************

void ext_main(void *r)
{
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("mda_stereo~", (method)mda_stereo_new, (method)mda_stereo_free, (long)sizeof(t_mda_stereo), 0L, A_GIMME, 0);
        
    class_addmethod(c, (method)mda_stereo_width, "width", A_GIMME, 0);
    class_addmethod(c, (method)mda_stereo_delay, "delay", A_GIMME, 0);
    class_addmethod(c, (method)mda_stereo_balance, "balance", A_GIMME, 0);
    class_addmethod(c, (method)mda_stereo_mod, "mod", A_GIMME, 0);
    class_addmethod(c, (method)mda_stereo_rate, "rate", A_GIMME, 0);
    class_addmethod(c, (method)mda_stereo_mode, "mode", A_GIMME, 0);
    class_addmethod(c, (method)mda_stereo_reset, "reset", A_GIMME, 0);
        
	class_addmethod(c, (method)mda_stereo_dsp64,    "dsp64",	A_CANT, 0);
	class_addmethod(c, (method)mda_stereo_assist,	"assist",	A_CANT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	mda_stereo_class = c;
}


void *mda_stereo_new(t_symbol *s, long argc, t_atom *argv)
{
	t_mda_stereo *x = (t_mda_stereo *)object_alloc(mda_stereo_class);

	if (x) {
        dsp_setup((t_pxobject *)x, 2);    // MSP inlets: arg is # of inlets and is REQUIRED!
                                
		// use 0 if you don't need inlets
                
		outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
        outlet_new(x, "signal");         // signal outlet (note "signal" rather than NULL)
    
        x->size = 4800;
        x->bufpos = 0;
        
        x->fParam1 = 0.78;
        x->fParam2 = 0.43;
        x->fParam3 = 0.50;
        x->fParam4 = 0.00;
        x->fParam5 = 0.50;
        x->mda_stereo_mode = gensym("comb");
                
        x->buffer = malloc(sizeof(double)*x->size);
        memset(x->buffer, 0, sizeof(double)*x->size);
        x->phi = 0;
        
        mda_stereo_recalculate(x);
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void mda_stereo_free(t_mda_stereo *x)
{
    dsp_free((t_pxobject *) x);
	
    free(x->buffer);
}


void mda_stereo_assist(t_mda_stereo *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) { //inlet
		sprintf(s, "I am inlet %ld", a);
	}
	else {	// outlet
		sprintf(s, "I am outlet %ld", a);
	}
}

void mda_stereo_width(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv)
{
    t_float f = atom_getfloat(argv);
    if(f < 0.0 || f > 1.0) {
        post("width should be between 0.0 and 1.0");
        return;
    }
    x->fParam1 = f;
    mda_stereo_recalculate(x);
}

void mda_stereo_delay(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv)
{
    t_float f = atom_getfloat(argv);
    if(f < 0.0 || f > 1.0) {
        post("delay should be between 0.0 and 1.0");
        return;
    }
    x->fParam2 = f;
    mda_stereo_recalculate(x);
}

void mda_stereo_balance(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv)
{
    t_float f = atom_getfloat(argv);
    if(f < 0.0 || f > 1.0) {
        post("balance should be between 0.0 and 1.0");
        return;
    }
    x->fParam3 = f;
    mda_stereo_recalculate(x);
}

void mda_stereo_mod(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv)
{
    t_float f = atom_getfloat(argv);
    if(f < 0.0 || f > 1.0) {
        post("mod should be between 0.0 and 1.0");
        return;
    }
    x->fParam4 = f;
    mda_stereo_recalculate(x);
}

void mda_stereo_rate(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv)
{
    t_float f = atom_getfloat(argv);
    if(f < 0.0 || f > 1.0){
        post("rate should be between 0.0 and 1.0");
        return;
    }
    x->fParam5 = f;
    mda_stereo_recalculate(x);
}

void mda_stereo_mode(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv)
{
    t_symbol *mda_stereo_mode = atom_getsym(argv);
    if(mda_stereo_mode == gensym("haas") || mda_stereo_mode == gensym("comb")){
        x->mda_stereo_mode = mda_stereo_mode;
        mda_stereo_recalculate(x);
        return;
    }
    post("mode should be 'haas' or 'comb'");
}

void mda_stereo_reset(t_mda_stereo *x, t_symbol *s, long argc, t_atom *argv)
{
    memset(x->buffer, 0, sizeof(double)*x->size);
}

// registers a function for the signal chain in Max
void mda_stereo_dsp64(t_mda_stereo *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the arguments passed are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: the symbol of the "dsp_add64" message we are sending
	// 3: a pointer to your object
	// 4: a pointer to your 64-bit perform method
	// 5: flags to alter how the signal chain handles your object -- just pass 0
	// 6: a generic pointer that you can use to pass any additional data to your perform method

	object_method(dsp64, gensym("dsp_add64"), x, mda_stereo_perform64, 0, NULL);
}

void mda_stereo_recalculate(t_mda_stereo *x)
{
    x->dphi= (float)(3.141 * pow (10.0,-2.0 + 3.0 * x->fParam5) / sys_getsr());
    x->mod= (float)(2100.0 * pow (x->fParam4, 2));
        
    if (x->mda_stereo_mode == gensym("haas"))
    {
        x->fli = 1.0f - (x->fParam1 * 0.75f);
        x->fld = 0.0f;
        x->fri = (float)(1.0 - x->fParam1);
        x->frd = (float)(1.0 - x->fri);
    }
    else
    {
        x->fli = ((1.0f - x->fParam1) * 0.5f) + 0.5f;
        x->fld = x->fParam1 * 0.5f;
        x->fri = x->fli;
        x->frd = -x->fld;
    }
    x->fdel = (float)(20.0 + 2080.0 * pow (x->fParam2, 2));
    if (x->fParam3>0.5)
    {
        x->fli *= (float)((1.0 - x->fParam3) * 2.0);
        x->fld *= (float)((1.0 - x->fParam3) * 2.0);
    }
    else
    {
        x->fri *= (2 * x->fParam3);
        x->frd *= (2 * x->fParam3);
    }
    x->fri *= (float)(0.5 + fabs (x->fParam1 - 0.5));
    x->frd *= (float)(0.5 + fabs (x->fParam1 - 0.5));
    x->fli *= (float)(0.5 + fabs (x->fParam1 - 0.5));
    x->fld *= (float)(0.5 + fabs (x->fParam1 - 0.5));
}

// this is the 64-bit perform method audio vectors
void mda_stereo_perform64(t_mda_stereo *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    
    long sampleFrames = sampleframes;
    
    double* in1 = ins[0];
    double* in2 = ins[1];
    double* out1 = outs[0];
    double* out2 = outs[1];

    float a, b, c, d;
    float li, ld, ri, rd, del, ph=x->phi, dph=x->dphi, mo=x->mod;
    long tmp, bp = x->bufpos;

    li = x->fli;
    ld = x->fld;
    ri = x->fri;
    rd = x->frd;
    del = x->fdel;

    --in1;
    --in2;
    --out1;
    --out2;
    if (mo>0.f) //modulated delay
    {
        while (--sampleFrames >= 0)
        {
            a = *++in1 + *++in2; //sum to mono

            *(x->buffer + bp) = a; //write
            tmp = (bp + (int)(del + fabs (mo * sin(ph)) ) ) % 4410;
            b = *(x->buffer + tmp);

            c = (a * li) - (b * ld); // output
            d = (a * ri) - (b * rd);

            bp = (bp - 1); if (bp < 0) bp = 4410; //buffer position

            ph = ph + dph;

            *++out1 = c;
            *++out2 = d;
        }
    }
    else
    {
        while (--sampleFrames >= 0)
        {
            a = *++in1 + *++in2; //sum to mono

            *(x->buffer + bp) = a; //write
            tmp = (bp + (int)(del) ) % 4410;
            b = *(x->buffer + tmp);

            c = (a * li) - (b * ld); // output
            d = (a * ri) - (b * rd);

            bp = (bp - 1); if (bp < 0) bp = 4410; //buffer position

            *++out1 = c;
            *++out2 = d;
        }
    }
    x->bufpos = bp;
    x->phi = (float)fmod(ph,6.2831853f);
}

