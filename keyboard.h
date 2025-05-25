#ifndef KEYBOARD_H
#define KEYBOARD_H

#ifndef __GNUC__
#include <conio.h>
#define INITKEYBOARD
#define DEINITKEYBOARD
#else
#define INITKEYBOARD ttycbreak()
#define DEINITKEYBOARD ttynorm()

/*   Keyboard support functions. January 7 1993 */

/* obligatory includes. */
#include <stdio.h>
#include <sys/types.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <termios.h>
#include <sys/ioctl.h>                                  
#ifndef FIONREAD
#define FIONREAD TIOCINQ
#endif
#include <signal.h>

static struct termios itio, tio;

static int ttyfd  = -1;
static int gottio = 0;

/* Catch quit signals so we dont end up with a CBREAK terminal when
   someone does a ^c during getch() */
void sig1 (int sig)
{
    tcsetattr (ttyfd, TCSANOW, &itio);
    signal (sig, SIG_DFL);
    kill (getpid(), sig);
}

/* initialize signals */
void setsigs()
{  
    signal (SIGINT,  sig1);
    signal (SIGTSTP, sig1);
    signal (SIGQUIT, sig1);
}

/* remove signals */
void rmsigs()
{
    signal (SIGCONT, SIG_DFL);
    signal (SIGINT,  SIG_DFL);
    signal (SIGTSTP, SIG_DFL);
    signal (SIGQUIT, SIG_DFL);
}

/* Put terminal in CBREAK mode */
int ttycbreak()
{
    if (ttyfd == -1)
	if ((ttyfd = open("/dev/tty", O_RDWR)) < 0)
	    return 0;
    if (tcgetattr(ttyfd, &tio) < 0)
	return 0;
    itio = tio;
    setsigs();
    gottio = 1;
    tio.c_cc[VMIN] = 0;
    tio.c_cc[VTIME] = 0;
    tio.c_lflag &= ~(ECHO|ICANON);
    tcsetattr (ttyfd, TCSANOW, &tio);
    return 1;
}

/* Restore terminal */
int ttynorm()
{
    gottio = 0;
    tcsetattr (ttyfd, TCSANOW, &itio);
    rmsigs();
    return 1;
}

/* Unix kbhit() */
int kbhit()
{
    long lf;
    if (ioctl(ttyfd, FIONREAD, &lf) == -1)
	return 0;
    return lf;
}

/* Unix getch() */
unsigned char getch()
{
    char c;
    read(ttyfd, &c, 1);
    return (unsigned char) c;
}

#endif
#endif
