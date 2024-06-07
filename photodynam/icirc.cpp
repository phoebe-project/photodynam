/*****************************************************************************/
/* icirc.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A library for computing intersections of circles.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2011; Pal, A. (apal@szofi.net)					     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "icirc.h"

int main(int argc,char *argv[])
{
 FILE	*fr,*fw;
 circle	*circles;
 int	ncircle;
 char	buff[256];
 arc	*arcs;
 int	narc;
 int	i;

 circles=NULL;
 ncircle=0;

 fr=stdin;
 while ( ! feof(fr) )
  {	double	x,y,r;
	if ( fgets(buff,256,fr)==NULL )
		break;
	if ( sscanf(buff,"%lg %lg %lg",&x,&y,&r)<3 )
		continue;
	circles=(circle *)realloc(circles,sizeof(circle)*(ncircle+1));
	(circles+ncircle)->x0=x;
	(circles+ncircle)->y0=y;
	(circles+ncircle)->r =r;
	ncircle++;
  }

 icirc_arclist_intersections(circles,ncircle,&arcs,&narc);

 fw=stdout;

 if ( 0 )
  {	for ( i=0 ; i<narc ; i++ )
	 {	arc	*a;
		int	j;
		a=&arcs[i];
		fprintf(fw,"%11g %11g %11g %11g %11g #",
			circles[a->cidx].x0,circles[a->cidx].y0,circles[a->cidx].r,
			a->phi0,a->dphi);
		for ( j=0 ; j<a->noidx && a->oidxs != NULL ; j++ )
	 	 {	fprintf(fw," %d",a->oidxs[j]);		}
		fprintf(fw,"\n");
	 }
  }
 if ( 1 )
  {	for ( i=0 ; i<narc ; i++ )
	 {	arc	*a;
		circle	*c;
		int	l,n;
		a=&arcs[i];
		n=100;
		for ( l=0 ; l<=n ; l++ )
		 {	double	x,y,p;
			c=&circles[a->cidx];
			p=a->phi0+a->dphi*(double)l/(double)n;
			x=c->x0+c->r*cos(p);
			y=c->y0+c->r*sin(p);
			fprintf(fw,"%d %11g %11g\n",(a->cidx==0&&a->oidxs==NULL?1:a->cidx&&a->oidxs&&a->noidx==1&&a->oidxs[0]==0?2:0),x,y);
		 }
		fprintf(fw,"\n");
	  }
  }

 return(0);
}

