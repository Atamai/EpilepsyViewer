//
//  Document.h
//  Epilepsy Viewer
//
//  Created by Yves Starreveld on 12-01-23.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@class BasicVTKView;

@interface Document : NSDocument
{
    IBOutlet BasicVTKView* mainVTKView;
}

@end
