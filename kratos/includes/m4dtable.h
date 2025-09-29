//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Marx Xu
//                   
//

#if !defined(KRATOS_M4DTABLE_H_INCLUDED )
#define  KRATOS_M4DTABLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "input_output/logger.h"
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/** 
 * @class M4dTable
 * @ingroup KratosCore
 * @brief This class represents the value of its variable depending to other variables.
 * @details M4dTable class stores the value of its second variable respect to the value of its first variable.
 * It also provides a double to double table with piecewise linear interpolator/extrapolator for getting intermediate values.
 * @author Marx Xu
 */
class M4dTable
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of M4dTable
    KRATOS_CLASS_POINTER_DEFINITION(M4dTable);

    /// Type definition of M4dTable
    typedef std::size_t RowType;
    typedef std::size_t ColumnType;
    typedef std::array<ColumnType, 4> ColumnListType;
    typedef std::vector<ColumnType> ColumnContainerType;
    
    typedef double ValueType;
    typedef std::array<ValueType, 4> ValueListType; 

    typedef std::vector<ValueType> ValueContainerType;
    typedef std::pair<std::string, ValueContainerType> RecordType;
    typedef std::array<RecordType, 4> RecordListType;

    ///@}
    ///@name Life Cycle
    ///@{

    M4dTable(): mDim(0), mSZ{0,0,0,0} {}

    M4dTable(RowType const& m): mDim(m), mSZ{0,0,0,0} {}

    virtual ~M4dTable() = default;

    ///@}
    ///@name Operators
    ///@{

    // This operator return the value of specific columns of argument variables. 
    ValueType& operator()(ColumnListType const& nn)
    {  
        switch(mDim)
        {
            case 1:

                return mData[nn[0]];
                break;

            case 2:

                return mData[mSZ[0]*nn[1] + nn[0]];
                break;

            case 3:

                return mData[mSZ[0]*(mSZ[1]*nn[2]+ nn[1]) + nn[0]];
                break;

            case 4:

                return mData[mSZ[0]*(mSZ[1]*(mSZ[2]*nn[3] + nn[2]) + nn[1]) + nn[0]];
                break;

            default:

        	return mData[0];
        }
    }

    // This operator return the value of specific columns of argument variables. 
    ValueType const& operator()(ColumnListType const& nn) const
    {  
        switch(mDim)
        {
            case 1:

                return mData[nn[0]];
                break;

            case 2:

                return mData[mSZ[0]*nn[1] + nn[0]];
                break;

            case 3:

                return mData[mSZ[0]*(mSZ[1]*nn[2]+ nn[1]) + nn[0]];
                break;

            case 4:

                return mData[mSZ[0]*(mSZ[1]*(mSZ[2]*nn[3] + nn[2]) + nn[1]) + nn[0]];
                break;

            default:

        	return mData[0];
        }
    }

    // This operator return the value of specific columns of argument variables. 
    ValueType& operator[](ColumnListType const& nn)
    {  
        switch(mDim)
        {
            case 1:

                return mData[nn[0]];
                break;

            case 2:

                return mData[mSZ[0]*nn[1] + nn[0]];
                break;

            case 3:

                return mData[mSZ[0]*(mSZ[1]*nn[2]+ nn[1]) + nn[0]];
                break;

            case 4:

                return mData[mSZ[0]*(mSZ[1]*(mSZ[2]*nn[3] + nn[2]) + nn[1]) + nn[0]];
                break;

            default:

        	return mData[0];
        }
    }

    // This operator return the value of specific columns of argument variables. 
    ValueType const& operator[](ColumnListType const& nn) const
    {  
        switch(mDim)
        {
            case 1:

                return mData[nn[0]];
                break;

            case 2:

                return mData[mSZ[0]*nn[1] + nn[0]];
                break;

            case 3:

                return mData[mSZ[0]*(mSZ[1]*nn[2]+ nn[1]) + nn[0]];
                break;

            case 4:

                return mData[mSZ[0]*(mSZ[1]*(mSZ[2]*nn[3] + nn[2]) + nn[1]) + nn[0]];
                break;

            default:

        	return mData[0];
        }
    }

    ///@}
    ///@name Operations
    ///@{

    // Get the nearest row according to the given value of argument 
    void GetNearestColumn(RowType const& m, ValueType const& x, ColumnType& nl, ColumnType& nr) const
    {
        // get the argument data object and size
        const ValueContainerType& arg = mArgument[m].second ;
        ColumnType size = arg.size();
      
        KRATOS_ERROR_IF(size == 0) << "Get row from empty argument table" << std::endl;

        // constant argument data
        if(size == 1)
        {
            nr = 0;
            nl = 0;
            return;
        }

        // check the x value
        if(x <= arg[0])
        {
            nr = 1;
            nl = 0;
            return;
        }

        for(ColumnType j = 1; j < size; j++)
        {
            if(x <= arg[j])
            {
                nr = j;
                nl = j - 1;
                return;
            }
        }

        nr = size - 1;
        nl = nr - 1;
        return;
    } 

    // Get the interpolate value according to the given value of argument 
    ValueType GetInterpolateValue(ValueListType const& x) const
    {
        ColumnListType nn, nl, nr;

        KRATOS_ERROR_IF(mDim == 0) << "Get interpolate value from empty m4dtable" << std::endl;

        // return the value on the recursive mode
        if(mDim > 0)
        {
            // get the nearest column information
            for(RowType m = 0; m < mDim; m++)
                GetNearestColumn(m, x[m], nl[m], nr[m]);

            return interpolate(0, x, nn, nl, nr);
        }
        else
        {
            return 0;
        }
    }

    // Get the derivative value according to the given value of argument 
    ValueType GetDerivativeValue(RowType const& m, ValueListType const& x) const
    {
        ColumnListType nn, nl, nr;
        ValueType xx1, xx2;
        ValueType f1, f2;

        KRATOS_ERROR_IF(mDim == 0) << "Get derivative value from empty m4dtable" << std::endl;

        // return the value on the recursive mode
        if(mDim > 0)
        {
            // copy the list of x
            ValueListType xx(x);

            // get the nearest column information
            for(RowType i = 0; i < mDim; i++)
                GetNearestColumn(i, x[i], nl[i], nr[i]);

            // check the nearest column
            if(nl[m] != nr[m])
            {
                xx1 = mArgument[m].second[nl[m]];
                xx[m] = xx1;
                f1 = interpolate(0, xx, nn, nl, nr);

                xx2 = mArgument[m].second[nr[m]];
                xx[m] = xx2;
                f2 = interpolate(0, xx, nn, nl, nr);
 
                return ((xx2 != xx1)? (f2 - f1) / (xx2 - xx1) : 0);
            }
            else
            {
                return 0;
            }
        }
        else
        {
            return 0;
        }
    }

    // Get the integrate value according to the given value of argument 
    ValueType GetIntegrationValue(ValueListType const& x1, ValueListType const& x2) const
    {
        ColumnListType nn, nl1, nr1, nl2, nr2;

        KRATOS_ERROR_IF(mDim == 0) << "Get integrate value from empty m4dtable" << std::endl;

        // return the value on the recursive mode
        if(mDim > 0 && IsNotEqual(x1, x2))
        {
            // get the nearest column information
            for(RowType m = 0; m < mDim; m++)
            {
                GetNearestColumn(m, x1[m], nl1[m], nr1[m]);
                GetNearestColumn(m, x2[m], nl2[m], nr2[m]);
            }

            return integrate(0, x1, x2, nn, nl1, nr1, nl2, nr2);
        }
        else
        {
            return 0;
        }
    }

    // compare two value list, if they are different, return True, otherwise return False 
    bool IsNotEqual(ValueListType const& x1, ValueListType const& x2) const
    {
        bool neq = false;

        for (unsigned int i = 0; i < mDim; i++)
        {
            if(x1[i] != x2[i]) {
                neq = true;
                break;
            }
        }              

        return neq;
    }

    // sort the data
    void Sort()
    {
        std::array<ColumnContainerType, 4> nn;
        ColumnType n1[3], n2[3];

        // sort the argument
        for(RowType m = 0; m < mDim; m++)
            value_sort(mArgument[m].second, nn[m]);

        ValueContainerType tmp(mData);

        // sort the data
        switch(mDim)
        {
            case 1:        

                for(ColumnType i = 0; i < mArgument[0].second.size(); i++)
                    tmp[i] = mData[nn[0][i]];

                mData = tmp;
                break;

            case 2:

                for(ColumnType i = 0; i < mArgument[1].second.size(); i++) {
                    n1[0] = mSZ[0]*i;
                    n2[0] = mSZ[0]*nn[1][i];

                    for(ColumnType j = 0; j < mArgument[0].second.size(); j++)
                        tmp[n1[0] + j] = mData[n2[0] + nn[0][j]];
                }

                mData = tmp;
                break;

            case 3:

                for(ColumnType i = 0; i < mArgument[2].second.size(); i++) {
                    n1[0] = mSZ[1]*i;
                    n2[0] = mSZ[1]*nn[2][i];

                    for(ColumnType j = 0; j < mArgument[1].second.size(); j++) {
                        n1[1] = mSZ[0]*(n1[0] + j);
                        n2[1] = mSZ[0]*(n2[0] + nn[1][j]);

                        for(ColumnType k = 0; k < mArgument[0].second.size(); k++)
                            tmp[n1[1] + k] = mData[n2[1] + nn[0][k]];
                    }
                }

                mData = tmp;
                break;

            case 4:

                for(ColumnType i = 0; i < mArgument[3].second.size(); i++) {
                    n1[0] = mSZ[2]*i;
                    n2[0] = mSZ[2]*nn[3][i];

                    for(ColumnType j = 0; j < mArgument[2].second.size(); j++) {
                        n1[1] = mSZ[1]*(n1[0] + j);
                        n2[1] = mSZ[1]*(n2[0] + nn[2][j]);

                        for(ColumnType k = 0; k < mArgument[1].second.size(); k++) {
                            n1[2] = mSZ[0]*(n1[1] + k);
                            n2[2] = mSZ[0]*(n2[1] + nn[1][k]);

                            for(ColumnType l = 0; l < mArgument[0].second.size(); l++)
                                tmp[n1[2] + l] = mData[n2[2] + nn[0][l]];
                        }
                    }
                }

                mData = tmp;
                break;

            default:

                break;
        }        
    }

    // put the argument row at the end.
    void PushBack(RowType const& m, ValueType const& x)
    {
        // check the dim
        if(m < mDim)
        {
            mArgument[m].second.push_back(x);
            mSZ[m] = mArgument[m].second.size();
        }
    }

    // put the data row at the end.
    void PushBack(ValueType const& x)
    {
        mData.push_back(x);
    }
    
    /**
     * @brief This method clears database
     */
    void Clear()
    {
        mData.clear();

        for(RowType i = 0; i < mDim; i++)
        {
            mArgument[i].first.clear();
            mArgument[i].second.clear();
        }
 
        mDim = 0;

        for(RowType j = 0; j < 4; j++)
            mSZ[j] = 0;
    }

    ///@}
    ///@name Access
    ///@{

    RecordListType& Argument()
    {
        return mArgument;
    }

    RecordListType const& Argument() const
    {
        return mArgument;
    }

    ValueContainerType& Data()
    {
        return mData;
    }

    ValueContainerType const& Data() const
    {
        return mData;
    }

    RowType& Dim()
    {
        return mDim;
    } 

    RowType const& Dim() const
    {
        return mDim;
    } 

    ColumnType const& Size(RowType const& m) const
    {
        if(m < mDim)
            return mSZ[m];

        return mSZ[0];
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "M4dTable";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << mNameOfY << "\t\t" << mDim << std::endl;

        for(RowType i = 0; i < mDim; i++)
            rOStream << mArgument[i].first << "\t\t" << mArgument[i].second.size() << std::endl;
    }

    const std::string& NameOfX(RowType const& m) const
    {
        // check the dim
        if(m < mDim)
            return mArgument[m].first;

        return mArgument[0].first;
    }

    const std::string& NameOfY() const
    {
        return mNameOfY;
    }

    void SetNameOfX(RowType const& m, const std::string& name)
    {
        if(m < mDim)
            mArgument[m].first = name;
    }

    void SetNameOfY(const std::string& name)
    {
        mNameOfY = name;
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    RowType mDim;
    ColumnType mSZ[4];
    ValueContainerType mData;
    RecordListType mArgument;
    std::string mNameOfY;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // sort the value container and return the new position information
    void value_sort(ValueContainerType& x, ColumnContainerType& n)
    {
        ColumnType local_size;
        ColumnType mn, tn;
        ValueType mv, tv;

        // adjust the size of n and init n
        local_size = x.size();
        n.resize(local_size);
        for(ColumnType i = 0; i < local_size; i++) 
            n[i] = i;

        for(ColumnType j = 0; j < local_size - 1; j++)
        {
            mn = j;
            mv = x[mn];
            
            // search the min value and position
            for(ColumnType k = j + 1; k < local_size; k++)
            {
                if(x[k] < mv)
                {
                    mn = k;
                    mv = x[mn];
                }
            }

            if(mn > j)
            {
                // swap the data                
                tv = x[j];
                x[j] = mv;
                x[mn] = tv;

                // update the position;
                tn = n[j];
                n[j] = n[mn];
                n[mn] = tn; 
            }
        }
    }

    // get the interpolation value with recursive mode
    ValueType interpolate(RowType const& m, ValueListType const& x, ColumnListType& nn, ColumnListType const& nl, ColumnListType const& nr) const
    {
        ValueType f1, f2;
        ValueType x1, x2;

        // check the dim
        if(m >= mDim)
        {
            return 0;
        }
        else if(m == (mDim - 1))
        {
            if(nl[m] == nr[m])
            {
                nn[m] = nl[m];
                return this->operator[](nn);  
            }
            else
            {
                nn[m] = nl[m];
                x1 = mArgument[m].second[nl[m]];
                f1= this->operator[](nn);

                nn[m] = nr[m];
                x2 = mArgument[m].second[nr[m]];
                f2 = this->operator[](nn);

                return (f2 - f1) * (x[m] - x1) / (x2 - x1) + f1;                
            }
        }
        else
        {
            if(nl[m] == nr[m])
            {
                nn[m] = nl[m];
                return this->interpolate(m+1, x, nn, nl, nr);  
            }
            else
            {
                nn[m] = nl[m];
                x1 = mArgument[m].second[nl[m]];
                f1= this->interpolate(m+1, x, nn, nl, nr);

                nn[m] = nr[m];
                x2 = mArgument[m].second[nr[m]];
                f2 = this->interpolate(m+1, x, nn, nl, nr);

                return (f2 - f1) * (x[m] - x1) / (x2 - x1) + f1;                
            }
        }
    }

    // get the integration value with recursive mode
    ValueType integrate(RowType const& m, ValueListType const& x1, ValueListType const& x2, ColumnListType& nn, ColumnListType const& nl1, ColumnListType const& nr1, ColumnListType const& nl2, ColumnListType const& nr2) const
    {
        ValueType f1, f2;
        ValueType xx1,xx2, xx;
        ValueType minx,maxx;
        ValueType sum = 0;
        ColumnType minn, maxn;
        ColumnType i, j;

        // get the integral scape on current dimension
        if(x1[m] < x2[m]) {
            minx = x1[m];
            maxx= x2[m];
        } else {
            minx = x2[m];
            maxx= x1[m];
        }

        // get the argument object
        ValueContainerType const& px = mArgument[m].second;

        // check the dim
        if(m >= mDim) {
            return 0;
        } else if(m == (mDim - 1)) {
            // check the integral scape
            if(minx == maxx) {
                if(nl1[m] == nr1[m]) {
                    nn[m] = nl1[m];                    
                    return this->operator[](nn);  
                } else {   
                    minn =  nl1[m];
                    xx1 = px[minn];
                    nn[m] = minn;
                    f1= this->operator[](nn);

                    maxn = nr1[m];
                    xx2 = px[maxn];
                    nn[m] = maxn;
                    f2 = this->operator[](nn);

                    return (f2 - f1) * (minx - xx1) / (xx2 - xx1) + f1;                
                }
            } else {
                if(mSZ[m] == 1) {
                    nn[m] = 0;
                    return (maxx - minx)*this->operator[](nn);
                } else {
                    minn = nl1[m] < nl2[m] ? nl1[m] : nl2[m];
                    maxn = nr1[m] < nr2[m] ? nr2[m] : nr1[m];
                    xx = minx;

                    for(i = minn; i < maxn; i++) {
                        j = i + 1;
                        xx2 = px[j];

                        if(xx < xx2) {
                            xx1 = px[i];
                            nn[m] = i;
                            f1 = this->operator[](nn);

                            nn[m] = j;
                            f2 = this->operator[](nn);

                            if(maxx < xx2) {
                                sum += ((f2- f1) * (xx + maxx) * 0.5 + xx2*f1 - xx1*f2) * (maxx - xx) / (xx2 - xx1);
                                return sum;
                            } else {
                                sum += ((f2- f1) * (xx + xx2) * 0.5 + xx2*f1 - xx1*f2) * (xx2 - xx) / (xx2 - xx1);
                                xx = xx2;
                            }
                        }
                    }

                    // check the rest part
                    if(xx < maxx) {
                        xx2 = px[i];
                        nn[m] = i;
                        f2 = this->operator[](nn);

                        xx1 = px[--i];
                        nn[m] = i;
                        f1 = this->operator[](nn);

                        sum += ((f2- f1) * (xx + maxx) * 0.5 + xx2*f1 - xx1*f2) * (maxx - xx) / (xx2 - xx1); 
                    }

                    return sum;
                }
            }
        } else {
            // check the integral scape
            if(minx == maxx) {
                if(nl1[m] == nr1[m]) {
                    nn[m] = nl1[m];    
                    return integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);
                } else {   
                    minn =  nl1[m];
                    xx1 = px[minn];
                    nn[m] = minn;
                    f1= integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);

                    maxn = nr1[m];
                    xx2 = px[maxn];
                    nn[m] = maxn;
                    f2 = integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);

                    return (f2 - f1) * (minx - xx1) / (xx2 - xx1) + f1;                
                }
            } else {
                if(mSZ[m] == 1) {
                    nn[m] = 0;
                    return (maxx - minx)*integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);
                } else {
                    minn = nl1[m] < nl2[m] ? nl1[m] : nl2[m];
                    maxn = nr1[m] < nr2[m] ? nr2[m] : nr1[m];
                    xx = minx;

                    for(i = minn; i < maxn; i++) {
                        j = i + 1;
                        xx2 = px[j];

                        if(xx < xx2) {
                            xx1 = px[i];
                            nn[m] = i;
                            f1 = integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);

                            nn[m] = j;
                            f2 = integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);

                            if(maxx < xx2) {
                                sum += ((f2- f1) * (xx + maxx) * 0.5 + xx2*f1 - xx1*f2) * (maxx - xx) / (xx2 - xx1);
                                return sum;
                            } else {
                                sum += ((f2- f1) * (xx + xx2) * 0.5 + xx2*f1 - xx1*f2) * (xx2 - xx) / (xx2 - xx1);
                                xx = xx2;
                            }
                        }
                    }

                    // check the rest part
                    if(xx < maxx) {
                        xx2 = px[i];
                        nn[m] = i;
                        f2 = integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);

                        xx1 = px[--i];
                        nn[m] = i;
                        f1 = integrate(m + 1, x1, x2, nn, nl1, nr1, nl2, nr2);

                        sum += ((f2- f1) * (xx + maxx) * 0.5 + xx2*f1 - xx1*f2) * (maxx - xx) / (xx2 - xx1);

                    }

                    return sum;
                }
            }           
        }

        return sum;
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {   
        rSerializer.save("name", mNameOfY);     
        rSerializer.save("dim", mDim);

        for(RowType i = 0; i < mDim; i++)
        {
            rSerializer.save("argument", mArgument[i].first);
            rSerializer.save("argu_size", mSZ[i]);

            for(auto j = mArgument[i].second.begin(); j != mArgument[i].second.end(); i++)
                rSerializer.save("argu_data", *j);
        }

        rSerializer.save("data_size", mData.size());

        for(auto k = mData.begin() ; k != mData.end() ; k++)
            rSerializer.save("data", *k);           
    }

    virtual void load(Serializer& rSerializer)
    {
        ColumnType local_size;

        rSerializer.load("name", mNameOfY);     
        rSerializer.load("dim", mDim);

        for(RowType i = 0; i < mDim; i++)
        {
            rSerializer.load("argument", mArgument[i].first);
            rSerializer.load("argu_size", mSZ[i]);
     
            mArgument[i].second.resize(mSZ[i]);

            for(auto j = mArgument[i].second.begin(); j != mArgument[i].second.end(); i++)
                rSerializer.load("argu_data", *j);
        }

        rSerializer.load("data_size", local_size);

        mData.resize(local_size);

        for(auto k = mData.begin() ; k != mData.end() ; k++)
            rSerializer.load("data", *k);     
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class M4dTable

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  M4dTable& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const M4dTable& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_M4DTABLE_H_INCLUDED  defined 