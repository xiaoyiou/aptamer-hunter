/*
 * chem.c
 *
 *  Created on: Jan 30, 2012
 *      Author: Yiou
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LINE_WIDTH 50
/**
 * need to avoid the dependencies on this number if they are trying to do
 * things on 64-char long reads.
 *
 */
#define READ_WIDTH 36
/*
 * Three categories
 */
#define GOOD 0
#define CANDIDATE 1
#define BAD 2
/*
 * The default size of the hash table. Later I might need to add more options
 * if the range of the hash function is enlarged (e.g. the sequence is 64 bits long)
 */
#define HASHP 100003
/**
 * Three report modes
 */
#define VERBOSE 2
#define CLUSTER 1
#define EXACT 0

/**
 * The node structure in the hash table
 */
typedef struct sequence {
	char *read;
	char *mid;
	struct sequence * next;
	int flag;
	int cnt;
	short ref;
} seq;
/**
 * this is just for the exact mode results printing
 */
typedef struct listnode {
	seq *read;
	struct listnode *next;
} lnode;

lnode goodHead;
lnode candHead;

/**
 * the global variables which needed to be initialized in the main
 * function later
 */
int headw = 10;
int tailw = 10;
int midw = 16;
int mismatch = 2;
int mode = EXACT;
int numOfLines = 0;
int color=0;

//Three tables, maybe I can get rid of the Bad table
seq *Good[HASHP];
seq *Cand[HASHP];
seq *Bad[HASHP];

//struct cluster *cls[HASHP]; //This one is not in use yet

/*Function heads*/
int distance(char *, char *, int, int);
int calIndex3(char *);
int classify(char *, char *, char *);

/**
 * This function might be used later to increase the efficiency
 *
 */
int calIndex3(char *str) {
	int sum = 0, i = 0;
	while (1) {
		switch (*(str + i)) {
		case 'A':
			break;
		case 'T':
			sum += 3;
			break;
		case 'C':
			sum += 2;
			break;
		case 'G':
			sum += 1;
			break;
		default:
			return sum;

		}
		i++;
	}
	return 0;
}

int classify(char *head, char *tail, char*read) {
	int hd, td;
	hd = distance(head, read, 0, headw);
	td = distance(tail, read, READ_WIDTH - tailw, tailw);
	if (hd > mismatch && td > mismatch)
		return BAD;
	else if (hd > mismatch || td > mismatch)
		return CANDIDATE;
	else
		return GOOD;
}

int calIndex2(char *str) {
	int i = 0, index = 0;
	while (1) {
		switch (*(str + i)) {
		case 'A':
		case 'T':
			index <<= 1;
			break;
		case 'C':
		case 'G':
			index = (index + 1) << 1;
			break;
		default:
			return index % HASHP;

		}
		i++;
	}
	return (index) % HASHP;
}

/**
 * this function will return the pointer to the target read in the table
 * However when the target is absent it will return a NULL pointer
 */
seq * find(seq *table[], char * target) {

	int index = calIndex2(target);
	seq * head = table[index];
	while (head != NULL) {
		if (strcmp(head->mid, target) == 0)
			return head;
		head = head->next;
	}
	return head;
}

void hashAdd(seq *table[], char *read) {

	seq* slot, *temp;
	int index, midw = READ_WIDTH - headw - tailw;
	char mid[(READ_WIDTH - headw - tailw) + 1];
	strncpy(mid, read+headw, midw);
	mid[midw] = '\0';
	index = calIndex2(mid);
	slot = find(table, mid);
	if (slot != NULL) {
		slot->cnt++;
		return;
	} else {
		temp = (seq *) malloc(sizeof(seq));
		temp->read = (char *) malloc(sizeof(char) * READ_WIDTH + 1);
		temp->mid = (char *) malloc(sizeof(char) * midw + 1);
		strcpy(temp->read, read);
		strncpy(temp->mid, mid, midw);
		temp->cnt = 1;
		temp->flag = 0;
		temp->next = NULL;
		slot = table[index];
		if (slot == NULL) {
			table[index] = temp;
		} else {
			while (slot->next != NULL)
				slot = slot->next;
			slot->next = temp;
		}
	}
}

/**
 * the l is the start point of the target fragment, w is the width of the
 * fragment
 * str1 is the fragment
 * str2 is the string need to be matched (must be longer than the fragment)
 */
int distance(char *str1, char *str2, int l, int w) {
	int i, dis = 0;
	for (i = 0; i < w; i++) {
		if (*(str1 + i) != *(str2 + i + l))
			dis++;
		if (dis > mismatch)
			break;
	}
	return dis;
}
/**
 * m is the criteria of # of acceptable mismatches.
 */
int readLines(char *filename, char * head, char *tail) {
	FILE* fp;
	int line = -1;
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot read %s\n", filename);
		return -1;
	}
	char *line_buf = (char *) malloc(sizeof(char) * LINE_WIDTH);
	while (fgets(line_buf, LINE_WIDTH, fp) != NULL) {
		line++;
		if (line % 4 == 1) {
			switch (classify(head, tail, line_buf)) {
			case GOOD:
				hashAdd(Good, line_buf);
				break;
			case CANDIDATE:
				hashAdd(Cand, line_buf);
				break;
			case BAD:
				hashAdd(Bad, line_buf);
				break;
			}
			numOfLines++;
		} else
			continue;
	}
	printf("The file has %d of sequences\n", numOfLines);
	return numOfLines;
}
void freeSeq(seq *p) {
	if (p == NULL)
		return;
	if (p->next != NULL) {
		freeSeq(p->next);
	}
	free(p->mid);
	free(p->read);
	free(p);
}
void freeLists(lnode* p) {
	if (p == NULL)
		return;
	if (p->next != NULL)
		freeLists(p->next);
	free(p);
}
void freeA() {
	int i;
	for (i = 0; i < HASHP; i++) {
		freeSeq(Good[i]);
		freeSeq(Cand[i]);
		freeSeq(Bad[i]);
	}
	freeLists(goodHead.next);
	freeLists(candHead.next);
}

int calMask(int origin, int pos) {
	int mask = 1;
	mask <<= pos;
	if ((origin >> pos) % 2 == 1) {
		origin = origin + mask - (mask << 1);
	} else
		origin = origin + mask;
	return origin;
}
void printM(char c) {
	if(color)
		printf("\e[31m%c\e[0m", c);
	else{
		printf("%c",c+32);
	}
}

int printBucket(seq *p, seq *m) {
	int k = 0, dis, i = 0;
	if (p == NULL) {
		return 0;
	} else {
		while (p != NULL) {

			dis = distance(m->mid, p->mid, 0, READ_WIDTH - tailw - headw);
			//printf("Target : %s\t This:%s mismatch:%d\n",m->mid,p->mid,dis);
			if (dis < 3 && dis > 0) {
				if (mode == VERBOSE) {
					printf("\tSequcne:"); //%s\tCount:%d\n", p->mid, p->cnt);
					for (i = 0; i < midw; i++) {
						if (p->mid[i] != m->mid[i])
							printM(p->mid[i]);
						else
							printf("%c", p->mid[i]);
					}
					printf(" Count:%d \n", p->cnt);
				}
				p->flag = 1;
				k += p->cnt;
			}
			p = p->next;

		}
	}
	return k;
}

void report(seq* a[], int threshold) {
	int i, k = 0, j = 0, jj = 0, total = 0;
	int index = 0, nindex = 0;
	seq *m;
	for (j = 0; j < HASHP; j++) {
		m = a[j];
		while (m) {
			total += m->cnt;
			m = m->next;
		}
	}
	/*
	 *  HERE I PRINT THE 2-nearest neighbors in the GOOD
	 */
	for (j = 0; j < HASHP; j++) {
		if (a[j] == NULL)
			continue;
		m = a[j];
		while (m) {
			if (m->cnt < threshold) {
				m = m->next;
				continue;
			}
			m->flag = 1;
			printf("==================\n\tSequnce:%s\tCount:%d\n", m->mid,
					m->cnt);
			//First 1 mismatch:
			index = calIndex2(m->mid);
			k += printBucket(a[j], m);
			for (i = 0; i < READ_WIDTH - headw - tailw; i++) {
				nindex = calMask(index, i);
				k += printBucket(a[nindex % HASHP], m);
			}
			//THE PROBLEM MIGHT BE HERE
			index = calIndex2(m->mid);
			for (i = 0; i < READ_WIDTH - headw - tailw - 1; i++) {
				for (jj = i + 1; jj < READ_WIDTH - headw - tailw; jj++) {
					nindex = calMask(calMask(index, i), jj);
					k += printBucket(a[nindex % HASHP], m);
				}
			}
			//Then 2 mismatches
			printf(
					"\tSize of this \"cluster\" %d\t %.2f%% of the total (total:%d)\n",
					k + m->cnt, (float) (k + m->cnt) / total * 100.0f, total);
			k = 0;
			m = m->next;
		}

	}
}
void insert(lnode *head, seq *target) {
	lnode *temp;
	temp = (lnode*) malloc(sizeof(lnode));
	temp->read = target;
	temp->next = NULL;
	while (head) {
		if (head->next == NULL) {
			head->next = temp;
			break;
		} else {
			if (head->next->read->cnt < target->cnt) {
				temp->next = head->next;
				head->next = temp;
				break;
			}
		}
		head = head->next;
	}

}
void printExact(int threshold) {
	int j = 0, i = 0;
	lnode *p1, *p2;
	seq *m, *goodTop, *candTop;
	//If I don't care about the cluster, only the exact occrences of each middle part!
	for (j = 0; j < HASHP; j++) {
		if (Good[j] == NULL)
			continue;
		m = Good[j];
		while (m) {
			if (m->cnt < threshold) {
				m = m->next;
				continue;
			}
			insert(&goodHead, m);
			m = m->next;
		}
	}

	for (j = 0; j < HASHP; j++) {
		if (Cand[j] == NULL)
			continue;
		m = Cand[j];
		while (m) {
			if (m->cnt < threshold) {
				m = m->next;
				continue;
			}
			insert(&candHead, m);
			m = m->next;
		}
	}

	p1 = goodHead.next;
	p2 = candHead.next;
	goodTop = p1->read;
	candTop = p2->read;
	printf("THE GOOD READS\t\t\t\t\t\tTHE CANDIDATE READS\n");
	while (p1 || p2) {
		if (p1) {

			//printf("Mid: %s count %d\t\t\t", p1->read->mid, p1->read->cnt);
			printf("Mid:");
			for (i = 0; i < midw; i++) {
				if (p1->read->mid[i] != goodTop->mid[i])
					printM(p1->read->mid[i]);
				else
					printf("%c", p1->read->mid[i]);
			}
			printf(" count %d\t\t\t", p1->read->cnt);
			p1 = p1->next;
		} else {
			printf("NOTHING\t\t\t\t\t\t");
		}
		if (p2) {
			printf("Mid:");
			for (i = 0; i < midw; i++) {
				if (p2->read->mid[i] != candTop->mid[i])
					printM(p2->read->mid[i]);
				else
					printf("%c", p2->read->mid[i]);
			}
			printf(" count %d\n", p2->read->cnt);
			p2 = p2->next;
		} else {
			printf("NOTHING\t\t\t\t\t\t\t");
		}
	}
}
void printUsage() {
	printf("\n\nSYNOPSIS\n\n");
	printf(
			"chem -i inputfile -h head -t tail [-C] [-T threshold] [-m mismatches] [-e|-v|-c]\n\n");
	printf("DESCRIPTION\n\n");
	printf(
			"-i inputfile ::\n\tis to define the path of the input file ILLUMIA format\n\n");
	printf(
			"[-T threshold] ::\n\tis optional to define the threshold the default value is 100\n\n");
	printf(
			"[-m mismatches] ::\n\tis to define the maximum number of mismatches acceptable\n\n");
	printf("[-C] ::\n\tcolor mode, which will only work in terminal not in file\n\n");
	printf(
			"[-e|-v|-c] ::\n\tis to define the report format -e is EXACT mode which means the exact counting,"
					"-c is CLUSTER mode which means the clustering mode without the details of each cluster -v is VERBOSE mode also"
					" clustering mode but it will print all the details\n");
}
void test() {
	int i = 0, j = 0, origin = 15;
	for (; i < 3; i++)
		for (j = i + 1; j < 4; j++)
			printf("Number %d \n", calMask(calMask(origin, i), j));
}
void print2b(int m, int bits) {
	bits--;
	if (bits > 0)
		print2b(m >> 1, bits);
	printf("%d", m % 2);
}
int main(int argc, char *argv[]) {
	int i, threshold = 100;
	//fprintf(stderr,"\033[;\033[0m");
	char *path = "/Users/xiaoyiou/s2", *head = "ACACGCGCATGC", *tail =
			"GCATGCGCC";

	 if (argc <= 1) {
	 printUsage();
	 return 0;
	 }
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 'i':
				path = argv[i + 1];
				break;
			case 'h':
				head = argv[i + 1];
				break;
			case 't':
				tail = argv[i + 1];
				break;
			case 'T':
				threshold = atoi(argv[i + 1]);
				break;
			case 'c':
				mode = CLUSTER;
				break;
			case 'e':
				mode = EXACT;
				break;
			case 'v':
				mode = VERBOSE;
				break;
			case 'm':
				mismatch = atoi(argv[i + 1]);
				break;
			case 'C':
				color=1;
				break;
			default:
				printUsage();
				return 0;
				break;
			}
		}
	}

	//TEST INPUT
	printf("The input File:%s, head %s, tail %s, threshold %d, mode %d, mismatches: %d\n",path,head,tail,threshold,mode,mismatch);

	headw = strlen(head);
	tailw = strlen(tail);
	midw = READ_WIDTH - headw - tailw;
	readLines(path, head, tail);
	if (mode > 0) {
		printf("Results of Good!\n");
		report(Good, threshold);
		printf("Results of Candidates\n");
		report(Good, threshold);
	} else {
		printExact(threshold);
	}
	freeA();

	return 0;

}
