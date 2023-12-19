#ifndef OPENMP_SETTINGS_H
#define OPENMP_SETTINGS_H

class OpenMPSettings {
private:
    int numThreads;

public:
    OpenMPSettings();

    int getNumThreads() const;

    void setNumThreads(int threads);
};

#endif