#ifndef TAGGEDVECTOR_H
#define TAGGEDVECTOR_H

#include <vector>
#include <queue>
#include <stdexcept>

class TaggedVector
{
    private:
        unsigned int numCurrEntries;
        unsigned int maxHistoricalTag;
        std::vector <int> v_tag;
        std::vector <int> v_rtag;
        std::queue <int> q_recycledTags;
        
    public:
        TaggedVector();
        
        int getNumEntries();
        
        int add();
		void remove(int oldTag);
        
        int getTagAtIndex(int pos);
		int getIndexOfTag(int tag);
        
        int getCurrMaxTag();
        
        void setTagAtPos(int pos, int newTag);
        bool checkTagExists(int tag);
};

#endif