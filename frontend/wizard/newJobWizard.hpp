/* 
 * File:   newJobWizard.hpp
 * Author: Petyr
 *
 * Created on 19. červen 2013, 17:01
 */

#ifndef _NEWJOBWIZARD_HPP
#define	_NEWJOBWIZARD_HPP

#include "ui_NewJobWizard.h"

class NewJobWizard : public QWizard {
    Q_OBJECT
public:
    NewJobWizard();
    virtual ~NewJobWizard();
private:
    Ui::newJobWizard ui;
};

#endif	/* _NEWJOBWIZARD_HPP */
