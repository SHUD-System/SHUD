#include <cstdio>
#include <cstring>

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Usage: %s <binary_file>\n", argv[0]);
        return 1;
    }
    
    const char* filename = argv[1];
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        printf("Error: Cannot open file %s\n", filename);
        return 1;
    }
    
    // Read header (1024 bytes)
    char header[1024];
    fread(header, sizeof(char), 1024, fp);
    
    // Read timestamp (double)
    double timestamp;
    fread(&timestamp, sizeof(double), 1, fp);
    
    // Read number of variables (double)
    double numVars;
    fread(&numVars, sizeof(double), 1, fp);
    
    fclose(fp);
    
    printf("File: %s\n", filename);
    printf("Header: %s\n", header);
    printf("Timestamp (modelBaseDate): %ld\n", (long)timestamp);
    printf("Number of variables: %d\n", (int)numVars);
    
    return 0;
}
